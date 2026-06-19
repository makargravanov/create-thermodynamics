use crate::chemistry::dynamic::DynamicChemistryRegistry;
use crate::chemistry::error::{ChemistryError, ChemistryResult};
use crate::chemistry::registry::ChemistryRegistry;

use super::blob::{
    decode_sections, encode_sections, NativeBlob, NativeBlobKind, NativeBlobLimits,
    NativeBlobSection,
};

pub const CATALOG_SNAPSHOT_MODEL_VERSION: &str = "create-thermodynamics:dynamic-catalog:1";

const SECTION_METADATA: u16 = 1;
const SECTION_DYNAMIC_SUBSTANCES: u16 = 2;

const CATALOG_METADATA_VERSION: u16 = 2;
const DYNAMIC_SUBSTANCES_VERSION: u16 = 1;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DynamicCatalogSnapshotMetadata {
    pub base_catalog_version: String,
    pub content_version: u64,
    pub static_substance_count: u64,
    pub static_reaction_count: u64,
    pub dynamic_substance_count: u64,
    pub dynamic_acid_base_spec_count: u64,
    pub dynamic_precipitation_spec_count: u64,
    pub dynamic_complex_spec_count: u64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DynamicSubstanceManifestEntry {
    pub substance_id: String,
    pub canonical_code: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DynamicCatalogSnapshotManifest {
    pub metadata: DynamicCatalogSnapshotMetadata,
    pub dynamic_substances: Vec<DynamicSubstanceManifestEntry>,
}

impl DynamicCatalogSnapshotManifest {
    pub fn from_registry(
        registry: &DynamicChemistryRegistry,
        base_catalog_version: impl Into<String>,
        content_version: u64,
    ) -> ChemistryResult<Self> {
        let base_catalog_version = base_catalog_version.into();
        validate_non_empty("base catalog version", &base_catalog_version)?;

        let static_registry = registry.static_registry();
        let dynamic_substances = registry
            .dynamic_substances()
            .map(|substance| DynamicSubstanceManifestEntry {
                substance_id: substance.id.as_str().to_string(),
                canonical_code: registry
                    .canonical_code_for(&substance.id)
                    .map(str::to_string),
            })
            .collect::<Vec<_>>();
        Ok(Self {
            metadata: DynamicCatalogSnapshotMetadata {
                base_catalog_version,
                content_version,
                static_substance_count: static_registry.substance_count() as u64,
                static_reaction_count: static_registry.reactions().count() as u64,
                dynamic_substance_count: dynamic_substances.len() as u64,
                dynamic_acid_base_spec_count: registry.dynamic_acid_base_spec_count() as u64,
                dynamic_precipitation_spec_count: registry.dynamic_precipitation_spec_count()
                    as u64,
                dynamic_complex_spec_count: registry.dynamic_complex_spec_count() as u64,
            },
            dynamic_substances,
        })
    }

    fn validate(&self) -> ChemistryResult<()> {
        validate_non_empty("base catalog version", &self.metadata.base_catalog_version)?;
        if self.metadata.dynamic_substance_count != self.dynamic_substances.len() as u64 {
            return Err(snapshot_error(format!(
                "dynamic substance manifest contains {} entries, metadata declares {}",
                self.dynamic_substances.len(),
                self.metadata.dynamic_substance_count
            )));
        }
        for entry in &self.dynamic_substances {
            validate_non_empty("dynamic substance id", &entry.substance_id)?;
            if let Some(canonical_code) = &entry.canonical_code {
                validate_non_empty("dynamic substance canonical code", canonical_code)?;
            }
        }
        Ok(())
    }
}

pub fn export_dynamic_catalog_checkpoint(
    registry: &DynamicChemistryRegistry,
    base_catalog_version: impl Into<String>,
    content_version: u64,
    limits: &NativeBlobLimits,
) -> ChemistryResult<Vec<u8>> {
    let manifest = DynamicCatalogSnapshotManifest::from_registry(
        registry,
        base_catalog_version,
        content_version,
    )?;
    manifest.validate()?;
    let payload = encode_manifest_sections(&manifest)?;
    NativeBlob::encode(
        NativeBlobKind::DynamicCatalogCheckpoint,
        CATALOG_SNAPSHOT_MODEL_VERSION,
        content_version,
        &payload,
        limits,
    )
}

pub fn read_dynamic_catalog_checkpoint_manifest(
    encoded: &[u8],
    expected_base_catalog_version: &str,
    limits: &NativeBlobLimits,
) -> ChemistryResult<DynamicCatalogSnapshotManifest> {
    validate_non_empty(
        "expected base catalog version",
        expected_base_catalog_version,
    )?;
    let (blob, payload) =
        NativeBlob::decode_expected(encoded, NativeBlobKind::DynamicCatalogCheckpoint, limits)?;
    if blob.model_version != CATALOG_SNAPSHOT_MODEL_VERSION {
        return Err(snapshot_error(format!(
            "dynamic catalog snapshot has model version '{}', expected '{}'",
            blob.model_version, CATALOG_SNAPSHOT_MODEL_VERSION
        )));
    }
    let manifest = decode_manifest_sections(&payload, limits)?;
    manifest.validate()?;
    if manifest.metadata.content_version != blob.content_version {
        return Err(snapshot_error(format!(
            "dynamic catalog metadata content version {} does not match blob content version {}",
            manifest.metadata.content_version, blob.content_version
        )));
    }
    if manifest.metadata.base_catalog_version != expected_base_catalog_version {
        return Err(snapshot_error(format!(
            "dynamic catalog snapshot is for base catalog '{}', expected '{}'",
            manifest.metadata.base_catalog_version, expected_base_catalog_version
        )));
    }
    Ok(manifest)
}

pub fn import_dynamic_catalog_checkpoint(
    base_registry: ChemistryRegistry,
    encoded: &[u8],
    expected_base_catalog_version: &str,
    limits: &NativeBlobLimits,
) -> ChemistryResult<DynamicChemistryRegistry> {
    let manifest =
        read_dynamic_catalog_checkpoint_manifest(encoded, expected_base_catalog_version, limits)?;
    let mut registry = DynamicChemistryRegistry::from_registry(base_registry)?;
    let static_registry = registry.static_registry();
    if manifest.metadata.static_substance_count != static_registry.substance_count() as u64 {
        return Err(snapshot_error(format!(
            "dynamic catalog snapshot declares {} static substances, current catalog has {}",
            manifest.metadata.static_substance_count,
            static_registry.substance_count()
        )));
    }
    if manifest.metadata.static_reaction_count != static_registry.reactions().count() as u64 {
        return Err(snapshot_error(format!(
            "dynamic catalog snapshot declares {} static reactions, current catalog has {}",
            manifest.metadata.static_reaction_count,
            static_registry.reactions().count()
        )));
    }
    for entry in &manifest.dynamic_substances {
        let canonical = entry.canonical_code.as_deref().ok_or_else(|| {
            snapshot_error(format!(
                "dynamic substance '{}' has no canonical code and cannot be restored",
                entry.substance_id
            ))
        })?;
        let restored_id = registry.resolve_frowns(canonical)?;
        if restored_id.as_str() != entry.substance_id {
            return Err(snapshot_error(format!(
                "dynamic substance canonical code '{}' restored '{}', expected '{}'",
                canonical, restored_id, entry.substance_id
            )));
        }
    }

    if manifest.metadata.dynamic_substance_count != registry.dynamic_substances().count() as u64 {
        return Err(snapshot_error(format!(
            "restored {} dynamic substances, snapshot declares {}",
            registry.dynamic_substances().count(),
            manifest.metadata.dynamic_substance_count
        )));
    }
    Ok(registry)
}

fn encode_manifest_sections(manifest: &DynamicCatalogSnapshotManifest) -> ChemistryResult<Vec<u8>> {
    encode_sections(&[
        NativeBlobSection {
            section_kind: SECTION_METADATA,
            payload: encode_metadata(&manifest.metadata)?,
        },
        NativeBlobSection {
            section_kind: SECTION_DYNAMIC_SUBSTANCES,
            payload: encode_dynamic_substances(&manifest.dynamic_substances)?,
        },
    ])
}

fn decode_manifest_sections(
    payload: &[u8],
    limits: &NativeBlobLimits,
) -> ChemistryResult<DynamicCatalogSnapshotManifest> {
    let mut metadata = None;
    let mut dynamic_substances = None;

    for section in decode_sections(payload, limits)? {
        match section.section_kind {
            SECTION_METADATA => assign_section(
                &mut metadata,
                decode_metadata(&section.payload)?,
                "metadata",
            )?,
            SECTION_DYNAMIC_SUBSTANCES => assign_section(
                &mut dynamic_substances,
                decode_dynamic_substances(&section.payload)?,
                "dynamic substances",
            )?,
            other => {
                return Err(snapshot_error(format!(
                    "unknown dynamic catalog snapshot section {other}"
                )));
            }
        }
    }

    Ok(DynamicCatalogSnapshotManifest {
        metadata: metadata
            .ok_or_else(|| snapshot_error("dynamic catalog snapshot is missing metadata"))?,
        dynamic_substances: dynamic_substances.ok_or_else(|| {
            snapshot_error("dynamic catalog snapshot is missing dynamic substance manifest")
        })?,
    })
}

fn assign_section<T>(target: &mut Option<T>, value: T, name: &str) -> ChemistryResult<()> {
    if target.is_some() {
        return Err(snapshot_error(format!(
            "dynamic catalog snapshot has duplicate {name} section"
        )));
    }
    *target = Some(value);
    Ok(())
}

fn encode_metadata(metadata: &DynamicCatalogSnapshotMetadata) -> ChemistryResult<Vec<u8>> {
    let mut output = Vec::new();
    write_u16(&mut output, CATALOG_METADATA_VERSION);
    write_string(&mut output, &metadata.base_catalog_version)?;
    write_u64(&mut output, metadata.content_version);
    write_u64(&mut output, metadata.static_substance_count);
    write_u64(&mut output, metadata.static_reaction_count);
    write_u64(&mut output, metadata.dynamic_substance_count);
    write_u64(&mut output, metadata.dynamic_acid_base_spec_count);
    write_u64(&mut output, metadata.dynamic_precipitation_spec_count);
    write_u64(&mut output, metadata.dynamic_complex_spec_count);
    Ok(output)
}

fn decode_metadata(payload: &[u8]) -> ChemistryResult<DynamicCatalogSnapshotMetadata> {
    let mut cursor = Cursor::new(payload);
    let version = cursor.read_u16()?;
    if version != CATALOG_METADATA_VERSION {
        return Err(snapshot_error(format!(
            "unsupported dynamic catalog metadata version {version}"
        )));
    }
    let metadata = DynamicCatalogSnapshotMetadata {
        base_catalog_version: cursor.read_string()?,
        content_version: cursor.read_u64()?,
        static_substance_count: cursor.read_u64()?,
        static_reaction_count: cursor.read_u64()?,
        dynamic_substance_count: cursor.read_u64()?,
        dynamic_acid_base_spec_count: cursor.read_u64()?,
        dynamic_precipitation_spec_count: cursor.read_u64()?,
        dynamic_complex_spec_count: cursor.read_u64()?,
    };
    cursor.finish()?;
    Ok(metadata)
}

fn encode_dynamic_substances(
    entries: &[DynamicSubstanceManifestEntry],
) -> ChemistryResult<Vec<u8>> {
    let mut output = Vec::new();
    write_u16(&mut output, DYNAMIC_SUBSTANCES_VERSION);
    write_u64(&mut output, entries.len() as u64);
    for entry in entries {
        write_string(&mut output, &entry.substance_id)?;
        write_optional_string(&mut output, entry.canonical_code.as_deref())?;
    }
    Ok(output)
}

fn decode_dynamic_substances(
    payload: &[u8],
) -> ChemistryResult<Vec<DynamicSubstanceManifestEntry>> {
    let mut cursor = Cursor::new(payload);
    let version = cursor.read_u16()?;
    if version != DYNAMIC_SUBSTANCES_VERSION {
        return Err(snapshot_error(format!(
            "unsupported dynamic substance manifest version {version}"
        )));
    }
    let count = cursor.read_len()?;
    let mut entries = Vec::with_capacity(count);
    for _ in 0..count {
        entries.push(DynamicSubstanceManifestEntry {
            substance_id: cursor.read_string()?,
            canonical_code: cursor.read_optional_string()?,
        });
    }
    cursor.finish()?;
    Ok(entries)
}

fn write_u16(output: &mut Vec<u8>, value: u16) {
    output.extend_from_slice(&value.to_le_bytes());
}

fn write_u64(output: &mut Vec<u8>, value: u64) {
    output.extend_from_slice(&value.to_le_bytes());
}

fn write_string(output: &mut Vec<u8>, value: &str) -> ChemistryResult<()> {
    validate_non_empty("string", value)?;
    write_u64(output, value.len() as u64);
    output.extend_from_slice(value.as_bytes());
    Ok(())
}

fn write_optional_string(output: &mut Vec<u8>, value: Option<&str>) -> ChemistryResult<()> {
    match value {
        Some(value) => {
            output.push(1);
            write_string(output, value)
        }
        None => {
            output.push(0);
            Ok(())
        }
    }
}

fn validate_non_empty(name: &str, value: &str) -> ChemistryResult<()> {
    if value.trim().is_empty() {
        return Err(snapshot_error(format!("{name} must not be empty")));
    }
    Ok(())
}

fn snapshot_error(reason: impl Into<String>) -> ChemistryError {
    ChemistryError::InvalidMixtureState(format!(
        "invalid dynamic catalog snapshot: {}",
        reason.into()
    ))
}

struct Cursor<'a> {
    bytes: &'a [u8],
    offset: usize,
}

impl<'a> Cursor<'a> {
    fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, offset: 0 }
    }

    fn read_u16(&mut self) -> ChemistryResult<u16> {
        let bytes = self.read_array::<2>()?;
        Ok(u16::from_le_bytes(bytes))
    }

    fn read_u64(&mut self) -> ChemistryResult<u64> {
        let bytes = self.read_array::<8>()?;
        Ok(u64::from_le_bytes(bytes))
    }

    fn read_len(&mut self) -> ChemistryResult<usize> {
        usize::try_from(self.read_u64()?)
            .map_err(|_| snapshot_error("length does not fit in usize"))
    }

    fn read_string(&mut self) -> ChemistryResult<String> {
        let len = self.read_len()?;
        let bytes = self.read_slice(len)?;
        let value = String::from_utf8(bytes.to_vec())
            .map_err(|error| snapshot_error(format!("invalid UTF-8 string: {error}")))?;
        validate_non_empty("string", &value)?;
        Ok(value)
    }

    fn read_optional_string(&mut self) -> ChemistryResult<Option<String>> {
        let tag = self.read_byte()?;
        match tag {
            0 => Ok(None),
            1 => Ok(Some(self.read_string()?)),
            other => Err(snapshot_error(format!(
                "invalid optional string marker {other}"
            ))),
        }
    }

    fn read_byte(&mut self) -> ChemistryResult<u8> {
        Ok(self.read_slice(1)?[0])
    }

    fn read_array<const N: usize>(&mut self) -> ChemistryResult<[u8; N]> {
        let slice = self.read_slice(N)?;
        let mut output = [0_u8; N];
        output.copy_from_slice(slice);
        Ok(output)
    }

    fn read_slice(&mut self, len: usize) -> ChemistryResult<&'a [u8]> {
        let end = self
            .offset
            .checked_add(len)
            .ok_or_else(|| snapshot_error("offset overflow"))?;
        if end > self.bytes.len() {
            return Err(snapshot_error(format!(
                "snapshot ended early at byte {}, needed {end}",
                self.bytes.len()
            )));
        }
        let slice = &self.bytes[self.offset..end];
        self.offset = end;
        Ok(slice)
    }

    fn finish(self) -> ChemistryResult<()> {
        if self.offset != self.bytes.len() {
            return Err(snapshot_error(format!(
                "snapshot has {} trailing bytes",
                self.bytes.len() - self.offset
            )));
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::dynamic::DynamicChemistryRegistry;

    fn limits() -> NativeBlobLimits {
        NativeBlobLimits::new(1024 * 1024, 8 * 1024 * 1024)
    }

    #[test]
    fn empty_dynamic_catalog_checkpoint_round_trips_manifest() {
        let registry = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let encoded =
            export_dynamic_catalog_checkpoint(&registry, "destroy-static:test", 1, &limits())
                .unwrap();

        let manifest =
            read_dynamic_catalog_checkpoint_manifest(&encoded, "destroy-static:test", &limits())
                .unwrap();

        assert_eq!(
            manifest.metadata.base_catalog_version,
            "destroy-static:test"
        );
        assert_eq!(manifest.metadata.content_version, 1);
        assert_eq!(manifest.metadata.dynamic_substance_count, 0);
        assert!(manifest.dynamic_substances.is_empty());
    }

    #[test]
    fn dynamic_catalog_checkpoint_records_dynamic_substances_without_reaction_cache() {
        let mut registry = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let dynamic_substance = registry.resolve_frowns("CCCCCCCC").unwrap();
        registry
            .generate_reactions_for(&dynamic_substance, 1)
            .unwrap();

        let encoded =
            export_dynamic_catalog_checkpoint(&registry, "destroy-static:test", 2, &limits())
                .unwrap();
        let manifest =
            read_dynamic_catalog_checkpoint_manifest(&encoded, "destroy-static:test", &limits())
                .unwrap();

        assert!(manifest.metadata.dynamic_substance_count > 0);
        assert_eq!(
            manifest.metadata.dynamic_substance_count,
            manifest.dynamic_substances.len() as u64
        );
        assert!(manifest
            .dynamic_substances
            .iter()
            .any(|entry| entry.substance_id == dynamic_substance.as_str()));
    }

    #[test]
    fn dynamic_catalog_checkpoint_rejects_wrong_base_catalog_version() {
        let registry = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let encoded =
            export_dynamic_catalog_checkpoint(&registry, "destroy-static:test", 1, &limits())
                .unwrap();

        let error =
            read_dynamic_catalog_checkpoint_manifest(&encoded, "destroy-static:other", &limits())
                .unwrap_err();

        assert!(error
            .to_string()
            .contains("expected 'destroy-static:other'"));
    }

    #[test]
    fn dynamic_catalog_checkpoint_import_restores_canonical_dynamic_substances() {
        let mut registry = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let dynamic_substance = registry.resolve_frowns("CCCCCCCC").unwrap();
        let encoded =
            export_dynamic_catalog_checkpoint(&registry, "destroy-static:test", 3, &limits())
                .unwrap();

        let restored = import_dynamic_catalog_checkpoint(
            crate::chemistry::destroy_registry_builder()
                .unwrap()
                .build()
                .unwrap(),
            &encoded,
            "destroy-static:test",
            &limits(),
        )
        .unwrap();

        assert!(restored.substance(&dynamic_substance).is_ok());
        assert_eq!(restored.dynamic_substances().count(), 1);
    }

    #[test]
    fn dynamic_catalog_checkpoint_import_discards_derived_reaction_cache() {
        let mut registry = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let dynamic_substance = registry.resolve_frowns("CCCCCCCC").unwrap();
        registry
            .generate_reactions_for(&dynamic_substance, 1)
            .unwrap();
        let encoded =
            export_dynamic_catalog_checkpoint(&registry, "destroy-static:test", 4, &limits())
                .unwrap();

        let restored = import_dynamic_catalog_checkpoint(
            crate::chemistry::destroy_registry_builder()
                .unwrap()
                .build()
                .unwrap(),
            &encoded,
            "destroy-static:test",
            &limits(),
        )
        .unwrap();

        assert!(restored.substance(&dynamic_substance).is_ok());
        assert_eq!(restored.dynamic_reactions().count(), 0);
    }

    #[test]
    fn dynamic_catalog_checkpoint_rejects_unknown_section() {
        let metadata = DynamicCatalogSnapshotMetadata {
            base_catalog_version: "destroy-static:test".to_string(),
            content_version: 1,
            static_substance_count: 0,
            static_reaction_count: 0,
            dynamic_substance_count: 0,
            dynamic_acid_base_spec_count: 0,
            dynamic_precipitation_spec_count: 0,
            dynamic_complex_spec_count: 0,
        };
        let payload = encode_sections(&[
            NativeBlobSection {
                section_kind: SECTION_METADATA,
                payload: encode_metadata(&metadata).unwrap(),
            },
            NativeBlobSection {
                section_kind: 999,
                payload: Vec::new(),
            },
            NativeBlobSection {
                section_kind: SECTION_DYNAMIC_SUBSTANCES,
                payload: encode_dynamic_substances(&[]).unwrap(),
            },
        ])
        .unwrap();
        let encoded = NativeBlob::encode(
            NativeBlobKind::DynamicCatalogCheckpoint,
            CATALOG_SNAPSHOT_MODEL_VERSION,
            1,
            &payload,
            &limits(),
        )
        .unwrap();

        let error =
            read_dynamic_catalog_checkpoint_manifest(&encoded, "destroy-static:test", &limits())
                .unwrap_err();

        assert!(error
            .to_string()
            .contains("unknown dynamic catalog snapshot section 999"));
    }
}
