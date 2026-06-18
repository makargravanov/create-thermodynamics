import { inflateSync, deflateSync } from "node:zlib";
import { mkdirSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { dirname } from "node:path";

const namespace = "create_thermodynamics";
const blockId = "reactor_chamber";
const resourceRoot = process.argv[2] ?? "modules/v1_21_1/v1_21_1-neoforge/src/generated/resources/assets/create_thermodynamics";
const blockbenchModelPath = "blockbench/reactor_chamber.bbmodel";
const connectedTemplatePath = "blockbench/reactor_chamber_connected_template.png";

const faceNames = ["north", "east", "south", "west", "up", "down"];

const blockbenchModel = JSON.parse(readFileSync(blockbenchModelPath, "utf8"));
const element = blockbenchModel.elements?.[0];
if (!element) {
  throw new Error(`${blockbenchModelPath} does not contain a block element`);
}

const baseTexture = blockbenchModel.textures?.find((texture) => texture.id === "0");
if (!baseTexture?.source?.startsWith("data:image/png;base64,")) {
  throw new Error(`${blockbenchModelPath} does not contain embedded base texture 0`);
}

const baseImage = decodePng(Buffer.from(baseTexture.source.slice("data:image/png;base64,".length), "base64"));
const connectedTemplate = decodePng(readFileSync(connectedTemplatePath));
if (baseImage.width !== 64 || baseImage.height !== 64) {
  throw new Error(`expected ${blockbenchModelPath} base texture to be 64x64`);
}
if (connectedTemplate.width !== 64 || connectedTemplate.height !== 64) {
  throw new Error(`expected ${connectedTemplatePath} to be 64x64`);
}

cleanupGeneratedChamberResources();
writeSplitBaseTextures(baseImage);
writeSplitConnectedTextures(connectedTemplate);
writeJson(`${resourceRoot}/models/block/${blockId}.json`, baseModel());
writeJson(`${resourceRoot}/models/item/${blockId}.json`, {
  parent: `${namespace}:block/${blockId}`,
});
writeJson(`${resourceRoot}/blockstates/${blockId}.json`, {
  variants: Object.fromEntries(
    chamberStates().map((state) => [
      chamberStateKey(state),
      { model: `${namespace}:block/${stateMask(state) === 0 ? blockId : `${blockId}_${stateMask(state).toString(16).padStart(2, "0")}`}` },
    ]),
  ),
});

for (const state of chamberStates().filter((state) => stateMask(state) !== 0)) {
  writeJson(
    `${resourceRoot}/models/block/${blockId}_${stateMask(state).toString(16).padStart(2, "0")}.json`,
    connectedModel(state),
  );
}

function cleanupGeneratedChamberResources() {
  rmSync(`${resourceRoot}/blockstates/${blockId}.json`, { force: true });
  rmSync(`${resourceRoot}/models/item/${blockId}.json`, { force: true });
  rmSync(`${resourceRoot}/textures/block/${blockId}.png`, { force: true });
  rmSync(`${resourceRoot}/textures/block/${blockId}_connected_template.png`, { force: true });
  rmSync(`${resourceRoot}/textures/block/${blockId}_connected.png`, { force: true });
  for (const face of faceNames) {
    rmSync(`${resourceRoot}/textures/block/${blockId}_${face}.png`, { force: true });
  }
  for (let mask = 0; mask < 16; mask += 1) {
    rmSync(`${resourceRoot}/textures/block/${blockId}_connected_${mask.toString(16)}.png`, { force: true });
    for (const face of faceNames) {
      rmSync(`${resourceRoot}/textures/block/${blockId}_connected_${face}_${mask.toString(16)}.png`, { force: true });
    }
  }
  for (let mask = 0; mask < 64; mask += 1) {
    rmSync(`${resourceRoot}/models/block/${blockId}_${mask.toString(16).padStart(2, "0")}.json`, { force: true });
  }
  rmSync(`${resourceRoot}/models/block/${blockId}.json`, { force: true });
}

function writeSplitBaseTextures(image) {
  for (const face of faceNames) {
    const data = element.faces[face];
    if (!data?.uv || data.uv.length !== 4) {
      throw new Error(`${blockbenchModelPath} face '${face}' does not contain a valid uv rectangle`);
    }
    writePng(`${resourceRoot}/textures/block/${blockId}_${face}.png`, cropUv(image, data.uv));
  }
}

function writeSplitConnectedTextures(image) {
  for (const face of faceNames) {
    const uv = element.faces[face]?.uv;
    if (!uv || uv.length !== 4) {
      throw new Error(`${blockbenchModelPath} face '${face}' does not contain a valid uv rectangle`);
    }
    const orientation = faceOrientation(uv);
    for (let mask = 0; mask < 16; mask += 1) {
      const sourceMask = orientMaskForSource(mask, orientation);
      const x = (sourceMask % 4) * 16;
      const y = Math.floor(sourceMask / 4) * 16;
      writePng(
        `${resourceRoot}/textures/block/${blockId}_connected_${face}_${mask.toString(16)}.png`,
        orientImage(cropRect(image, x, y, 16, 16), orientation),
      );
    }
  }
}

function baseModel() {
  return {
    textures: Object.fromEntries([
      ["particle", `${namespace}:block/${blockId}_north`],
      ...faceNames.map((face) => [face, `${namespace}:block/${blockId}_${face}`]),
    ]),
    elements: [
      {
        from: element.from,
        to: element.to,
        faces: Object.fromEntries(
          faceNames.map((face) => [
            face,
            {
              uv: [0, 0, 16, 16],
              texture: `#${face}`,
              cullface: face,
            },
          ]),
        ),
      },
    ],
  };
}

function connectedModel(state) {
  return {
    textures: {
      particle: `${namespace}:block/${blockId}_north`,
      north: connectedTexture("north", faceMask(state.up, state.east, state.down, state.west)),
      east: connectedTexture("east", faceMask(state.up, state.south, state.down, state.north)),
      south: connectedTexture("south", faceMask(state.up, state.west, state.down, state.east)),
      west: connectedTexture("west", faceMask(state.up, state.north, state.down, state.south)),
      up: connectedTexture("up", faceMask(state.north, state.east, state.south, state.west)),
      down: connectedTexture("down", faceMask(state.south, state.east, state.north, state.west)),
    },
    elements: [
      {
        from: element.from,
        to: element.to,
        faces: Object.fromEntries(
          faceNames.map((face) => [
            face,
            {
              uv: [0, 0, 16, 16],
              texture: `#${face}`,
              cullface: face,
            },
          ]),
        ),
      },
    ],
  };
}

function connectedTexture(face, mask) {
  return `${namespace}:block/${blockId}_connected_${face}_${mask.toString(16)}`;
}

function faceMask(up, right, down, left) {
  return (up ? 1 : 0) | (right ? 2 : 0) | (down ? 4 : 0) | (left ? 8 : 0);
}

function faceOrientation(uv) {
  return {
    flipX: uv[2] < uv[0],
    flipY: uv[3] < uv[1],
  };
}

function orientMaskForSource(mask, orientation) {
  let sourceMask = mask;
  if (orientation.flipX) {
    sourceMask = swapBits(sourceMask, 2, 8);
  }
  if (orientation.flipY) {
    sourceMask = swapBits(sourceMask, 1, 4);
  }
  return sourceMask;
}

function swapBits(mask, a, b) {
  const hasA = (mask & a) !== 0;
  const hasB = (mask & b) !== 0;
  mask &= ~a;
  mask &= ~b;
  if (hasA) mask |= b;
  if (hasB) mask |= a;
  return mask;
}

function orientImage(image, orientation) {
  if (!orientation.flipX && !orientation.flipY) {
    return image;
  }
  const output = emptyImage(image.width, image.height);
  for (let y = 0; y < image.height; y += 1) {
    for (let x = 0; x < image.width; x += 1) {
      const sourceX = orientation.flipX ? image.width - 1 - x : x;
      const sourceY = orientation.flipY ? image.height - 1 - y : y;
      copyPixel(image, output, sourceX, sourceY, x, y);
    }
  }
  return output;
}

function chamberStates() {
  const states = [];
  for (let mask = 0; mask < 64; mask += 1) {
    states.push({
      north: Boolean(mask & 1),
      east: Boolean(mask & 2),
      south: Boolean(mask & 4),
      west: Boolean(mask & 8),
      up: Boolean(mask & 16),
      down: Boolean(mask & 32),
    });
  }
  return states;
}

function stateMask(state) {
  return (
    (state.north ? 1 : 0) |
    (state.east ? 2 : 0) |
    (state.south ? 4 : 0) |
    (state.west ? 8 : 0) |
    (state.up ? 16 : 0) |
    (state.down ? 32 : 0)
  );
}

function chamberStateKey(state) {
  return [
    `north=${state.north}`,
    `east=${state.east}`,
    `south=${state.south}`,
    `west=${state.west}`,
    `up=${state.up}`,
    `down=${state.down}`,
  ].join(",");
}

function cropUv(image, uv) {
  const [u1, v1, u2, v2] = uv;
  const output = emptyImage(16, 16);
  for (let y = 0; y < 16; y += 1) {
    for (let x = 0; x < 16; x += 1) {
      const sourceX = Math.floor(u1 + ((u2 - u1) * (x + 0.5)) / 16);
      const sourceY = Math.floor(v1 + ((v2 - v1) * (y + 0.5)) / 16);
      copyPixel(image, output, clamp(sourceX, 0, image.width - 1), clamp(sourceY, 0, image.height - 1), x, y);
    }
  }
  return output;
}

function cropRect(image, x, y, width, height) {
  const output = emptyImage(width, height);
  for (let yy = 0; yy < height; yy += 1) {
    for (let xx = 0; xx < width; xx += 1) {
      copyPixel(image, output, x + xx, y + yy, xx, yy);
    }
  }
  return output;
}

function copyPixel(source, target, sourceX, sourceY, targetX, targetY) {
  const sourceIndex = (sourceY * source.width + sourceX) * 4;
  const targetIndex = (targetY * target.width + targetX) * 4;
  target.data[targetIndex] = source.data[sourceIndex];
  target.data[targetIndex + 1] = source.data[sourceIndex + 1];
  target.data[targetIndex + 2] = source.data[sourceIndex + 2];
  target.data[targetIndex + 3] = source.data[sourceIndex + 3];
}

function emptyImage(width, height) {
  return { width, height, data: new Uint8Array(width * height * 4) };
}

function decodePng(buffer) {
  const signature = "89504e470d0a1a0a";
  if (buffer.subarray(0, 8).toString("hex") !== signature) {
    throw new Error("invalid PNG signature");
  }

  let offset = 8;
  let width = 0;
  let height = 0;
  let colorType = 0;
  const idatChunks = [];
  while (offset < buffer.length) {
    const length = buffer.readUInt32BE(offset);
    const type = buffer.subarray(offset + 4, offset + 8).toString("ascii");
    const data = buffer.subarray(offset + 8, offset + 8 + length);
    offset += 12 + length;
    if (type === "IHDR") {
      width = data.readUInt32BE(0);
      height = data.readUInt32BE(4);
      const bitDepth = data[8];
      colorType = data[9];
      if (bitDepth !== 8 || colorType !== 6) {
        throw new Error(`unsupported PNG format: bitDepth=${bitDepth}, colorType=${colorType}`);
      }
    } else if (type === "IDAT") {
      idatChunks.push(data);
    } else if (type === "IEND") {
      break;
    }
  }

  if (!width || !height || colorType !== 6) {
    throw new Error("PNG does not contain a supported IHDR");
  }

  const bytesPerPixel = 4;
  const stride = width * bytesPerPixel;
  const inflated = inflateSync(Buffer.concat(idatChunks));
  const data = new Uint8Array(width * height * bytesPerPixel);
  let inputOffset = 0;
  let outputOffset = 0;
  const previous = new Uint8Array(stride);

  for (let y = 0; y < height; y += 1) {
    const filter = inflated[inputOffset++];
    const row = new Uint8Array(stride);
    row.set(inflated.subarray(inputOffset, inputOffset + stride));
    inputOffset += stride;
    unfilterRow(filter, row, previous, bytesPerPixel);
    data.set(row, outputOffset);
    previous.set(row);
    outputOffset += stride;
  }
  return { width, height, data };
}

function unfilterRow(filter, row, previous, bytesPerPixel) {
  for (let i = 0; i < row.length; i += 1) {
    const left = i >= bytesPerPixel ? row[i - bytesPerPixel] : 0;
    const up = previous[i];
    const upLeft = i >= bytesPerPixel ? previous[i - bytesPerPixel] : 0;
    switch (filter) {
      case 0:
        break;
      case 1:
        row[i] = (row[i] + left) & 0xff;
        break;
      case 2:
        row[i] = (row[i] + up) & 0xff;
        break;
      case 3:
        row[i] = (row[i] + Math.floor((left + up) / 2)) & 0xff;
        break;
      case 4:
        row[i] = (row[i] + paeth(left, up, upLeft)) & 0xff;
        break;
      default:
        throw new Error(`unsupported PNG filter: ${filter}`);
    }
  }
}

function paeth(left, up, upLeft) {
  const p = left + up - upLeft;
  const pa = Math.abs(p - left);
  const pb = Math.abs(p - up);
  const pc = Math.abs(p - upLeft);
  if (pa <= pb && pa <= pc) return left;
  if (pb <= pc) return up;
  return upLeft;
}

function encodePng(image) {
  const rows = new Uint8Array((image.width * 4 + 1) * image.height);
  for (let y = 0; y < image.height; y += 1) {
    rows[y * (image.width * 4 + 1)] = 0;
    rows.set(image.data.subarray(y * image.width * 4, (y + 1) * image.width * 4), y * (image.width * 4 + 1) + 1);
  }
  const ihdr = Buffer.alloc(13);
  ihdr.writeUInt32BE(image.width, 0);
  ihdr.writeUInt32BE(image.height, 4);
  ihdr[8] = 8;
  ihdr[9] = 6;
  ihdr[10] = 0;
  ihdr[11] = 0;
  ihdr[12] = 0;
  return Buffer.concat([
    Buffer.from([0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a]),
    chunk("IHDR", ihdr),
    chunk("IDAT", deflateSync(rows)),
    chunk("IEND", Buffer.alloc(0)),
  ]);
}

function chunk(type, data) {
  const typeBuffer = Buffer.from(type, "ascii");
  const length = Buffer.alloc(4);
  length.writeUInt32BE(data.length, 0);
  const crcBuffer = Buffer.alloc(4);
  crcBuffer.writeUInt32BE(crc(Buffer.concat([typeBuffer, data])), 0);
  return Buffer.concat([length, typeBuffer, data, crcBuffer]);
}

function crc(buffer) {
  let value = 0xffffffff;
  for (const byte of buffer) {
    value ^= byte;
    for (let bit = 0; bit < 8; bit += 1) {
      value = (value >>> 1) ^ (0xedb88320 & -(value & 1));
    }
  }
  return (value ^ 0xffffffff) >>> 0;
}

function writePng(path, image) {
  writeBytes(path, encodePng(image));
}

function writeJson(path, value) {
  writeText(path, `${JSON.stringify(value, null, 2)}\n`);
}

function writeText(path, text) {
  writeBytes(path, text);
}

function writeBytes(path, bytes) {
  mkdirSync(dirname(path), { recursive: true });
  writeFileSync(path, bytes);
}

function clamp(value, min, max) {
  return Math.max(min, Math.min(max, value));
}
