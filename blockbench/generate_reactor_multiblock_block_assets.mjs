import { deflateSync } from "node:zlib";
import { mkdirSync, writeFileSync } from "node:fs";
import { dirname } from "node:path";

const namespace = "create_thermodynamics";
const resourceRoot = "modules/v1_21_1/v1_21_1-neoforge/src/main/resources/assets/create_thermodynamics";

const blocks = [
  {
    id: "reactor_controller",
    displayName: "Reactor Controller",
    accent: rgb(183, 139, 62),
    marks: "controller",
  },
  {
    id: "reactor_item_input_port",
    displayName: "Reactor Item Input Port",
    accent: rgb(88, 143, 189),
    marks: "item_input",
  },
  {
    id: "reactor_item_output_port",
    displayName: "Reactor Item Output Port",
    accent: rgb(76, 175, 120),
    marks: "item_output",
  },
  {
    id: "reactor_fluid_input_port",
    displayName: "Reactor Fluid Input Port",
    accent: rgb(72, 164, 184),
    marks: "fluid_input",
  },
  {
    id: "reactor_fluid_output_port",
    displayName: "Reactor Fluid Output Port",
    accent: rgb(180, 91, 102),
    marks: "fluid_output",
  },
];

for (const block of blocks) {
  writePng(`${resourceRoot}/textures/block/${block.id}.png`, texture(block));
  writeText(`${resourceRoot}/models/block/${block.id}.json`, `${JSON.stringify({
    parent: "minecraft:block/cube_all",
    textures: {
      all: `${namespace}:block/${block.id}`,
    },
  }, null, 2)}\n`);
  writeText(`${resourceRoot}/models/item/${block.id}.json`, `${JSON.stringify({
    parent: `${namespace}:block/${block.id}`,
  }, null, 2)}\n`);
  writeText(`${resourceRoot}/blockstates/${block.id}.json`, `${JSON.stringify({
    variants: {
      "": {
        model: `${namespace}:block/${block.id}`,
      },
    },
  }, null, 2)}\n`);
  writeText(`blockbench/${block.id}.bbmodel`, `${JSON.stringify(blockbenchModel(block), null, 2)}\n`);
}

function texture(block) {
  const width = 16;
  const height = 16;
  const pixels = new Uint8Array(width * height * 4);
  const c = {
    frame: rgb(19, 23, 25),
    shadow: rgb(34, 39, 41),
    body: rgb(48, 54, 56),
    panel: rgb(73, 82, 85),
    panelLight: rgb(100, 112, 116),
    panelDark: rgb(38, 44, 46),
    black: rgb(11, 13, 14),
    white: rgb(205, 214, 210),
    accent: block.accent,
  };

  fill(pixels, width, 0, 0, 16, 16, c.body);
  rect(pixels, width, 0, 0, 16, 1, c.frame);
  rect(pixels, width, 0, 15, 16, 1, c.frame);
  rect(pixels, width, 0, 0, 1, 16, c.frame);
  rect(pixels, width, 15, 0, 1, 16, c.frame);
  rect(pixels, width, 2, 2, 12, 12, c.shadow);
  rect(pixels, width, 3, 3, 10, 10, c.panel);
  rect(pixels, width, 4, 4, 8, 1, c.panelLight);
  rect(pixels, width, 4, 11, 8, 1, c.panelDark);
  rect(pixels, width, 4, 4, 1, 8, c.panelLight);
  rect(pixels, width, 11, 4, 1, 8, c.panelDark);

  switch (block.marks) {
    case "controller":
      rect(pixels, width, 6, 5, 4, 1, c.accent);
      rect(pixels, width, 5, 6, 1, 4, c.accent);
      rect(pixels, width, 10, 6, 1, 4, c.accent);
      rect(pixels, width, 6, 10, 4, 1, c.accent);
      rect(pixels, width, 7, 7, 2, 2, c.black);
      point(pixels, width, 8, 8, c.white);
      break;
    case "item_input":
      arrow(pixels, width, c.accent, true);
      rect(pixels, width, 5, 5, 6, 1, c.black);
      rect(pixels, width, 5, 10, 6, 1, c.black);
      break;
    case "item_output":
      arrow(pixels, width, c.accent, false);
      rect(pixels, width, 5, 5, 6, 1, c.black);
      rect(pixels, width, 5, 10, 6, 1, c.black);
      break;
    case "fluid_input":
      drop(pixels, width, c.accent);
      arrow(pixels, width, c.white, true, 1);
      break;
    case "fluid_output":
      drop(pixels, width, c.accent);
      arrow(pixels, width, c.white, false, 1);
      break;
  }

  return encodePng(width, height, pixels);
}

function arrow(pixels, width, color, inward, yOffset = 0) {
  if (inward) {
    rect(pixels, width, 4, 8 + yOffset, 6, 1, color);
    rect(pixels, width, 8, 6 + yOffset, 1, 5, color);
    point(pixels, width, 9, 7 + yOffset, color);
    point(pixels, width, 10, 8 + yOffset, color);
    point(pixels, width, 9, 9 + yOffset, color);
  } else {
    rect(pixels, width, 6, 8 + yOffset, 6, 1, color);
    rect(pixels, width, 7, 6 + yOffset, 1, 5, color);
    point(pixels, width, 6, 7 + yOffset, color);
    point(pixels, width, 5, 8 + yOffset, color);
    point(pixels, width, 6, 9 + yOffset, color);
  }
}

function drop(pixels, width, color) {
  rect(pixels, width, 7, 4, 2, 1, color);
  rect(pixels, width, 6, 5, 4, 1, color);
  rect(pixels, width, 5, 6, 6, 3, color);
  rect(pixels, width, 6, 9, 4, 2, color);
  rect(pixels, width, 7, 11, 2, 1, color);
}

function blockbenchModel(block) {
  const textureUuid = uuidFor(`${block.id}-texture`);
  const elementUuid = uuidFor(`${block.id}-element`);
  return {
    meta: {
      format_version: "5.0",
      model_format: "java_block",
      box_uv: true,
    },
    name: block.displayName,
    model_identifier: `${namespace}:${block.id}`,
    visible_box: [1, 1, 1],
    variable_placeholders: "",
    variable_placeholder_buttons: [],
    timeline_setups: [],
    resolution: {
      width: 16,
      height: 16,
    },
    elements: [
      {
        name: block.id,
        box_uv: true,
        rescale: false,
        locked: false,
        from: [0, 0, 0],
        to: [16, 16, 16],
        autouv: 0,
        color: 0,
        origin: [8, 8, 8],
        faces: Object.fromEntries(["north", "east", "south", "west", "up", "down"].map((face) => [
          face,
          {
            uv: [0, 0, 16, 16],
            texture: 0,
          },
        ])),
        type: "cube",
        uuid: elementUuid,
      },
    ],
    outliner: [elementUuid],
    textures: [
      {
        path: `../modules/v1_21_1/v1_21_1-neoforge/src/main/resources/assets/${namespace}/textures/block/${block.id}.png`,
        name: block.id,
        folder: "block",
        namespace,
        id: "0",
        particle: true,
        render_mode: "default",
        render_sides: "auto",
        frame_time: 1,
        frame_order_type: "loop",
        frame_order: "",
        frame_interpolate: false,
        visible: true,
        mode: "link",
        saved: true,
        uuid: textureUuid,
        relative_path: `../modules/v1_21_1/v1_21_1-neoforge/src/main/resources/assets/${namespace}/textures/block/${block.id}.png`,
        source: "",
      },
    ],
    groups: [],
    animations: [],
    display: {},
  };
}

function writePng(path, bytes) {
  mkdirSync(dirname(path), { recursive: true });
  writeFileSync(path, bytes);
}

function writeText(path, text) {
  mkdirSync(dirname(path), { recursive: true });
  writeFileSync(path, text);
}

function fill(pixels, width, x, y, w, h, color) {
  rect(pixels, width, x, y, w, h, color);
}

function rect(pixels, width, x, y, w, h, color) {
  for (let yy = y; yy < y + h; yy += 1) {
    for (let xx = x; xx < x + w; xx += 1) {
      point(pixels, width, xx, yy, color);
    }
  }
}

function point(pixels, width, x, y, color) {
  if (x < 0 || y < 0 || x >= width || y >= 16) return;
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

function uuidFor(input) {
  const bytes = new Uint8Array(16);
  for (let i = 0; i < input.length; i += 1) {
    bytes[i % 16] = (bytes[i % 16] * 31 + input.charCodeAt(i)) & 0xff;
  }
  bytes[6] = (bytes[6] & 0x0f) | 0x40;
  bytes[8] = (bytes[8] & 0x3f) | 0x80;
  const hex = Array.from(bytes, (byte) => byte.toString(16).padStart(2, "0")).join("");
  return `${hex.slice(0, 8)}-${hex.slice(8, 12)}-${hex.slice(12, 16)}-${hex.slice(16, 20)}-${hex.slice(20)}`;
}
