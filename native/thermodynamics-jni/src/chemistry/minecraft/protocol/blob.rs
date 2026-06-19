use sha2::{Digest, Sha256};

use crate::chemistry::error::{ChemistryError, ChemistryResult};

const MAGIC: [u8; 8] = *b"CTNBLB1\0";
const FORMAT_VERSION: u16 = 1;
const HEADER_LEN: usize = 8 + 2 + 2 + 2 + 8 + 8 + 8 + 32;
const DEFAULT_ZSTD_LEVEL: i32 = 3;
const MAX_MODEL_VERSION_BYTES: usize = 256;
const MAX_SECTION_COUNT: usize = 65_536;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u16)]
pub enum NativeBlobKind {
    DynamicCatalogCheckpoint = 1,
    DynamicCatalogDelta = 2,
    ReactorSnapshot = 3,
    ReactionSpaceSnapshot = 4,
    PlannerSnapshot = 5,
}

impl NativeBlobKind {
    fn from_wire(value: u16) -> ChemistryResult<Self> {
        match value {
            1 => Ok(Self::DynamicCatalogCheckpoint),
            2 => Ok(Self::DynamicCatalogDelta),
            3 => Ok(Self::ReactorSnapshot),
            4 => Ok(Self::ReactionSpaceSnapshot),
            5 => Ok(Self::PlannerSnapshot),
            _ => Err(blob_error(format!("unknown native blob kind {value}"))),
        }
    }

    fn to_wire(self) -> u16 {
        self as u16
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NativeBlobLimits {
    pub max_compressed_bytes: usize,
    pub max_uncompressed_bytes: usize,
}

impl NativeBlobLimits {
    pub const fn new(max_compressed_bytes: usize, max_uncompressed_bytes: usize) -> Self {
        Self {
            max_compressed_bytes,
            max_uncompressed_bytes,
        }
    }
}

impl Default for NativeBlobLimits {
    fn default() -> Self {
        Self::new(64 * 1024 * 1024, 512 * 1024 * 1024)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NativeBlob {
    pub kind: NativeBlobKind,
    pub model_version: String,
    pub content_version: u64,
    pub content_hash: [u8; 32],
    pub uncompressed_len: usize,
    pub compressed_payload: Vec<u8>,
}

impl NativeBlob {
    pub fn encode(
        kind: NativeBlobKind,
        model_version: impl Into<String>,
        content_version: u64,
        payload: &[u8],
        limits: &NativeBlobLimits,
    ) -> ChemistryResult<Vec<u8>> {
        let model_version = model_version.into();
        validate_model_version(&model_version)?;
        validate_uncompressed_len(payload.len(), limits)?;

        let compressed_payload = zstd::bulk::compress(payload, DEFAULT_ZSTD_LEVEL)
            .map_err(|error| blob_error(format!("failed to compress native blob: {error}")))?;
        validate_compressed_len(compressed_payload.len(), limits)?;

        let content_hash = hash(payload);
        let model_version_bytes = model_version.as_bytes();
        let mut output =
            Vec::with_capacity(HEADER_LEN + model_version_bytes.len() + compressed_payload.len());
        output.extend_from_slice(&MAGIC);
        write_u16(&mut output, FORMAT_VERSION);
        write_u16(&mut output, kind.to_wire());
        write_u16(&mut output, model_version_bytes.len() as u16);
        write_u64(&mut output, content_version);
        write_u64(&mut output, payload.len() as u64);
        write_u64(&mut output, compressed_payload.len() as u64);
        output.extend_from_slice(&content_hash);
        output.extend_from_slice(model_version_bytes);
        output.extend_from_slice(&compressed_payload);
        Ok(output)
    }

    pub fn read(encoded: &[u8], limits: &NativeBlobLimits) -> ChemistryResult<Self> {
        let mut cursor = Cursor::new(encoded);
        let magic = cursor.read_array::<8>()?;
        if magic != MAGIC {
            return Err(blob_error("native blob has invalid magic"));
        }

        let format_version = cursor.read_u16()?;
        if format_version != FORMAT_VERSION {
            return Err(blob_error(format!(
                "unsupported native blob format version {format_version}"
            )));
        }

        let kind = NativeBlobKind::from_wire(cursor.read_u16()?)?;
        let model_version_len = cursor.read_u16()? as usize;
        if model_version_len > MAX_MODEL_VERSION_BYTES {
            return Err(blob_error(format!(
                "native blob model version is too long: {model_version_len} bytes"
            )));
        }
        let content_version = cursor.read_u64()?;
        let uncompressed_len = cursor.read_u64()? as usize;
        let compressed_len = cursor.read_u64()? as usize;
        validate_uncompressed_len(uncompressed_len, limits)?;
        validate_compressed_len(compressed_len, limits)?;
        let content_hash = cursor.read_array::<32>()?;
        let model_version = cursor.read_string(model_version_len)?;
        let compressed_payload = cursor.read_vec(compressed_len)?;
        cursor.finish()?;

        Ok(Self {
            kind,
            model_version,
            content_version,
            content_hash,
            uncompressed_len,
            compressed_payload,
        })
    }

    pub fn decode_payload(&self, limits: &NativeBlobLimits) -> ChemistryResult<Vec<u8>> {
        validate_uncompressed_len(self.uncompressed_len, limits)?;
        validate_compressed_len(self.compressed_payload.len(), limits)?;
        let payload = zstd::bulk::decompress(&self.compressed_payload, self.uncompressed_len)
            .map_err(|error| blob_error(format!("failed to decompress native blob: {error}")))?;
        if payload.len() != self.uncompressed_len {
            return Err(blob_error(format!(
                "native blob decompressed to {} bytes, expected {}",
                payload.len(),
                self.uncompressed_len
            )));
        }
        let actual_hash = hash(&payload);
        if actual_hash != self.content_hash {
            return Err(blob_error("native blob content hash mismatch"));
        }
        Ok(payload)
    }

    pub fn decode_expected(
        encoded: &[u8],
        expected_kind: NativeBlobKind,
        limits: &NativeBlobLimits,
    ) -> ChemistryResult<(Self, Vec<u8>)> {
        let blob = Self::read(encoded, limits)?;
        if blob.kind != expected_kind {
            return Err(blob_error(format!(
                "native blob has kind {:?}, expected {:?}",
                blob.kind, expected_kind
            )));
        }
        let payload = blob.decode_payload(limits)?;
        Ok((blob, payload))
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NativeBlobSection {
    pub section_kind: u16,
    pub payload: Vec<u8>,
}

pub fn encode_sections(sections: &[NativeBlobSection]) -> ChemistryResult<Vec<u8>> {
    if sections.len() > MAX_SECTION_COUNT {
        return Err(blob_error(format!(
            "native blob has too many sections: {}",
            sections.len()
        )));
    }
    let mut output = Vec::new();
    write_u32(&mut output, sections.len() as u32);
    for section in sections {
        if section.section_kind == 0 {
            return Err(blob_error("native blob section kind must not be zero"));
        }
        write_u16(&mut output, section.section_kind);
        write_u64(&mut output, section.payload.len() as u64);
        output.extend_from_slice(&section.payload);
    }
    Ok(output)
}

pub fn decode_sections(
    payload: &[u8],
    limits: &NativeBlobLimits,
) -> ChemistryResult<Vec<NativeBlobSection>> {
    let mut cursor = Cursor::new(payload);
    let section_count = cursor.read_u32()? as usize;
    if section_count > MAX_SECTION_COUNT {
        return Err(blob_error(format!(
            "native blob declares too many sections: {section_count}"
        )));
    }
    let mut sections = Vec::with_capacity(section_count);
    for _ in 0..section_count {
        let section_kind = cursor.read_u16()?;
        if section_kind == 0 {
            return Err(blob_error("native blob section kind must not be zero"));
        }
        let len = cursor.read_u64()? as usize;
        validate_uncompressed_len(len, limits)?;
        let payload = cursor.read_vec(len)?;
        sections.push(NativeBlobSection {
            section_kind,
            payload,
        });
    }
    cursor.finish()?;
    Ok(sections)
}

fn validate_model_version(model_version: &str) -> ChemistryResult<()> {
    if model_version.is_empty() {
        return Err(blob_error("native blob model version must not be empty"));
    }
    if model_version.len() > MAX_MODEL_VERSION_BYTES {
        return Err(blob_error(format!(
            "native blob model version is too long: {} bytes",
            model_version.len()
        )));
    }
    Ok(())
}

fn validate_uncompressed_len(len: usize, limits: &NativeBlobLimits) -> ChemistryResult<()> {
    if len > limits.max_uncompressed_bytes {
        return Err(blob_error(format!(
            "native blob uncompressed size {len} exceeds limit {}",
            limits.max_uncompressed_bytes
        )));
    }
    Ok(())
}

fn validate_compressed_len(len: usize, limits: &NativeBlobLimits) -> ChemistryResult<()> {
    if len > limits.max_compressed_bytes {
        return Err(blob_error(format!(
            "native blob compressed size {len} exceeds limit {}",
            limits.max_compressed_bytes
        )));
    }
    Ok(())
}

fn hash(payload: &[u8]) -> [u8; 32] {
    let mut hasher = Sha256::new();
    hasher.update(payload);
    hasher.finalize().into()
}

fn write_u16(output: &mut Vec<u8>, value: u16) {
    output.extend_from_slice(&value.to_le_bytes());
}

fn write_u32(output: &mut Vec<u8>, value: u32) {
    output.extend_from_slice(&value.to_le_bytes());
}

fn write_u64(output: &mut Vec<u8>, value: u64) {
    output.extend_from_slice(&value.to_le_bytes());
}

fn blob_error(reason: impl Into<String>) -> ChemistryError {
    ChemistryError::InvalidMixtureState(format!("invalid native blob: {}", reason.into()))
}

struct Cursor<'a> {
    bytes: &'a [u8],
    offset: usize,
}

impl<'a> Cursor<'a> {
    fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, offset: 0 }
    }

    fn read_array<const N: usize>(&mut self) -> ChemistryResult<[u8; N]> {
        let slice = self.read_slice(N)?;
        let mut output = [0_u8; N];
        output.copy_from_slice(slice);
        Ok(output)
    }

    fn read_u16(&mut self) -> ChemistryResult<u16> {
        Ok(u16::from_le_bytes(self.read_array::<2>()?))
    }

    fn read_u32(&mut self) -> ChemistryResult<u32> {
        Ok(u32::from_le_bytes(self.read_array::<4>()?))
    }

    fn read_u64(&mut self) -> ChemistryResult<u64> {
        Ok(u64::from_le_bytes(self.read_array::<8>()?))
    }

    fn read_string(&mut self, len: usize) -> ChemistryResult<String> {
        let bytes = self.read_slice(len)?;
        String::from_utf8(bytes.to_vec())
            .map_err(|error| blob_error(format!("invalid UTF-8 in model version: {error}")))
    }

    fn read_vec(&mut self, len: usize) -> ChemistryResult<Vec<u8>> {
        Ok(self.read_slice(len)?.to_vec())
    }

    fn read_slice(&mut self, len: usize) -> ChemistryResult<&'a [u8]> {
        let end = self
            .offset
            .checked_add(len)
            .ok_or_else(|| blob_error("native blob offset overflow"))?;
        if end > self.bytes.len() {
            return Err(blob_error(format!(
                "native blob ended early at byte {}, needed {end}",
                self.bytes.len()
            )));
        }
        let slice = &self.bytes[self.offset..end];
        self.offset = end;
        Ok(slice)
    }

    fn finish(self) -> ChemistryResult<()> {
        if self.offset != self.bytes.len() {
            return Err(blob_error(format!(
                "native blob has {} trailing bytes",
                self.bytes.len() - self.offset
            )));
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn small_limits() -> NativeBlobLimits {
        NativeBlobLimits::new(1024, 4096)
    }

    #[test]
    fn native_blob_round_trips_payload() {
        let payload = b"reactor snapshot payload";
        let encoded = NativeBlob::encode(
            NativeBlobKind::ReactorSnapshot,
            "test-model",
            7,
            payload,
            &small_limits(),
        )
        .unwrap();

        let (blob, decoded) =
            NativeBlob::decode_expected(&encoded, NativeBlobKind::ReactorSnapshot, &small_limits())
                .unwrap();

        assert_eq!(blob.kind, NativeBlobKind::ReactorSnapshot);
        assert_eq!(blob.model_version, "test-model");
        assert_eq!(blob.content_version, 7);
        assert_eq!(decoded, payload);
    }

    #[test]
    fn native_blob_rejects_wrong_kind() {
        let encoded = NativeBlob::encode(
            NativeBlobKind::DynamicCatalogCheckpoint,
            "test-model",
            1,
            b"catalog",
            &small_limits(),
        )
        .unwrap();

        let error =
            NativeBlob::decode_expected(&encoded, NativeBlobKind::ReactorSnapshot, &small_limits())
                .unwrap_err();

        assert!(error.to_string().contains("expected ReactorSnapshot"));
    }

    #[test]
    fn native_blob_rejects_modified_hash() {
        let mut encoded = NativeBlob::encode(
            NativeBlobKind::ReactorSnapshot,
            "test-model",
            1,
            b"payload",
            &small_limits(),
        )
        .unwrap();
        let hash_offset = 8 + 2 + 2 + 2 + 8 + 8 + 8;
        encoded[hash_offset] ^= 0x01;

        let blob = NativeBlob::read(&encoded, &small_limits()).unwrap();
        let error = blob.decode_payload(&small_limits()).unwrap_err();

        assert!(error.to_string().contains("hash mismatch"));
    }

    #[test]
    fn native_blob_rejects_oversized_payload_before_compression() {
        let limits = NativeBlobLimits::new(1024, 4);
        let error = NativeBlob::encode(
            NativeBlobKind::ReactorSnapshot,
            "test-model",
            1,
            b"too large",
            &limits,
        )
        .unwrap_err();

        assert!(error.to_string().contains("uncompressed size"));
    }

    #[test]
    fn native_blob_sections_round_trip() {
        let sections = vec![
            NativeBlobSection {
                section_kind: 1,
                payload: b"substances".to_vec(),
            },
            NativeBlobSection {
                section_kind: 2,
                payload: b"reactions".to_vec(),
            },
        ];
        let payload = encode_sections(&sections).unwrap();
        let encoded = NativeBlob::encode(
            NativeBlobKind::DynamicCatalogCheckpoint,
            "test-model",
            3,
            &payload,
            &small_limits(),
        )
        .unwrap();

        let (_, decoded_payload) = NativeBlob::decode_expected(
            &encoded,
            NativeBlobKind::DynamicCatalogCheckpoint,
            &small_limits(),
        )
        .unwrap();
        let decoded_sections = decode_sections(&decoded_payload, &small_limits()).unwrap();

        assert_eq!(decoded_sections, sections);
    }
}
