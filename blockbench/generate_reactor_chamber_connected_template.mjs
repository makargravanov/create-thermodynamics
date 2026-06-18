import { deflateSync } from "node:zlib";
import { writeFileSync } from "node:fs";

const width = 64;
const height = 64;
const tile = 16;

const colors = {
  body: rgb(47, 52, 54),
  bodyShadow: rgb(35, 40, 42),
  frame: rgb(21, 25, 27),
  frameInner: rgb(34, 40, 42),
  panel: rgb(73, 82, 85),
  panelLight: rgb(99, 112, 115),
  panelDark: rgb(37, 43, 45),
  brass: rgb(168, 135, 72),
};

const pixels = new Uint8Array(width * height * 4);
fillRect(0, 0, width, height, rgb(0, 0, 0, 0));

for (let mask = 0; mask < 16; mask += 1) {
  const x = (mask % 4) * tile;
  const y = Math.floor(mask / 4) * tile;
  drawTile(x, y, {
    up: (mask & 1) !== 0,
    right: (mask & 2) !== 0,
    down: (mask & 4) !== 0,
    left: (mask & 8) !== 0,
  }, mask);
}

writeOutput(
  "blockbench/reactor_chamber_connected_template.png",
  "blockbench/reactor_chamber_connected_template_fixed.png",
  encodePng(width, height, pixels),
);
writeOutput(
  "blockbench/reactor_chamber_connected_template.json",
  "blockbench/reactor_chamber_connected_template_fixed.json",
  `${JSON.stringify({
    tile_size: 16,
    atlas_size: [64, 64],
    mask_bits: {
      up: 1,
      right: 2,
      down: 4,
      left: 8,
    },
    tiles: Array.from({ length: 16 }, (_, mask) => ({
      mask,
      x: (mask % 4) * tile,
      y: Math.floor(mask / 4) * tile,
      connected: {
        up: (mask & 1) !== 0,
        right: (mask & 2) !== 0,
        down: (mask & 4) !== 0,
        left: (mask & 8) !== 0,
      },
    })),
  }, null, 2)}\n`,
);

function writeOutput(path, busyFallbackPath, content) {
  try {
    writeFileSync(path, content);
  } catch (error) {
    if (error?.code !== "EBUSY") {
      throw error;
    }
    writeFileSync(busyFallbackPath, content);
    console.warn(`${path} is locked; wrote ${busyFallbackPath} instead`);
  }
}

function drawTile(x, y, connected, seed) {
  fillRect(x, y, tile, tile, colors.body);

  if (!connected.up) fillRect(x, y, tile, 1, colors.frame);
  if (!connected.down) fillRect(x, y + tile - 1, tile, 1, colors.frame);
  if (!connected.left) fillRect(x, y, 1, tile, colors.frame);
  if (!connected.right) fillRect(x + tile - 1, y, 1, tile, colors.frame);

  if (!connected.up) fillRect(x + (connected.left ? 0 : 1), y + 1, tile - (connected.left ? 0 : 1) - (connected.right ? 0 : 1), 2, colors.bodyShadow);
  if (!connected.down) fillRect(x + (connected.left ? 0 : 1), y + tile - 3, tile - (connected.left ? 0 : 1) - (connected.right ? 0 : 1), 2, colors.bodyShadow);
  if (!connected.left) fillRect(x + 1, y + (connected.up ? 0 : 1), 2, tile - (connected.up ? 0 : 1) - (connected.down ? 0 : 1), colors.bodyShadow);
  if (!connected.right) fillRect(x + tile - 3, y + (connected.up ? 0 : 1), 2, tile - (connected.up ? 0 : 1) - (connected.down ? 0 : 1), colors.bodyShadow);

  const panelLeft = connected.left ? 0 : 3;
  const panelRight = connected.right ? tile - 1 : tile - 4;
  const panelTop = connected.up ? 0 : 3;
  const panelBottom = connected.down ? tile - 1 : tile - 4;
  fillRect(x + panelLeft, y + panelTop, panelRight - panelLeft + 1, panelBottom - panelTop + 1, colors.panel);

  if (!connected.left) fillRect(x + 2, y + panelTop, 1, panelBottom - panelTop + 1, colors.frameInner);
  if (!connected.right) fillRect(x + tile - 3, y + panelTop, 1, panelBottom - panelTop + 1, colors.frameInner);
  if (!connected.up) fillRect(x + panelLeft, y + 2, panelRight - panelLeft + 1, 1, colors.frameInner);
  if (!connected.down) fillRect(x + panelLeft, y + tile - 3, panelRight - panelLeft + 1, 1, colors.frameInner);

  if (!connected.left && !connected.up) setPixel(x + 1, y + 1, colors.brass);
  if (!connected.right && !connected.up) setPixel(x + tile - 2, y + 1, colors.brass);
  if (!connected.left && !connected.down) setPixel(x + 1, y + tile - 2, colors.brass);
  if (!connected.right && !connected.down) setPixel(x + tile - 2, y + tile - 2, colors.brass);

  const lightPoints = [
    [5, 4], [8, 4], [10, 7], [6, 11],
  ];
  const darkPoints = [
    [4, 8], [11, 5], [12, 10],
  ];
  for (const [px, py] of lightPoints) {
    const tx = x + ((px + seed) % tile);
    const ty = y + py;
    if (tx >= x + panelLeft && tx <= x + panelRight && ty >= y + panelTop && ty <= y + panelBottom) {
      setPixel(tx, ty, colors.panelLight);
    }
  }
  for (const [px, py] of darkPoints) {
    const tx = x + px;
    const ty = y + ((py + seed) % tile);
    if (tx >= x + panelLeft && tx <= x + panelRight && ty >= y + panelTop && ty <= y + panelBottom) {
      setPixel(tx, ty, colors.panelDark);
    }
  }
}

function fillRect(x, y, w, h, color) {
  for (let yy = y; yy < y + h; yy += 1) {
    for (let xx = x; xx < x + w; xx += 1) {
      setPixel(xx, yy, color);
    }
  }
}

function setPixel(x, y, color) {
  if (x < 0 || y < 0 || x >= width || y >= height) return;
  const index = (y * width + x) * 4;
  pixels[index] = color[0];
  pixels[index + 1] = color[1];
  pixels[index + 2] = color[2];
  pixels[index + 3] = color[3];
}

function rgb(r, g, b, a = 255) {
  return [r, g, b, a];
}

function encodePng(w, h, rgba) {
  const rows = new Uint8Array((w * 4 + 1) * h);
  for (let y = 0; y < h; y += 1) {
    rows[y * (w * 4 + 1)] = 0;
    rows.set(rgba.subarray(y * w * 4, (y + 1) * w * 4), y * (w * 4 + 1) + 1);
  }
  const ihdr = Buffer.alloc(13);
  ihdr.writeUInt32BE(w, 0);
  ihdr.writeUInt32BE(h, 4);
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
