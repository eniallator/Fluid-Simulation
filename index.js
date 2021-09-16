const canvas = document.getElementById("canvas");
const ctx = canvas.getContext("2d");

const mouse = new Mouse(canvas);
const paramConfig = new ParamConfig("./config.json", $("#cfg-outer")[0]);
paramConfig.addCopyToClipboardHandler("#share-btn");

window.onresize = (evt) => {
  canvas.width = $("#canvas").width();
  canvas.height = $("#canvas").height();
};
window.onresize();

ctx.fillStyle = "black";
ctx.strokeStyle = "white";

const fluid = new Fluid(256, 0.01, 0, 0);
let oldMousePos = mouse.relativePos.copy();

function clamp(val, min, max) {
  return Math.min(Math.max(val, min), max);
}

function run() {
  if (mouse.down) {
    fluid.addDensity(
      Math.floor(mouse.relativePos.x * Fluid.N),
      Math.floor(mouse.relativePos.y * Fluid.N),
      10000
    );
    const diff = mouse.relativePos.copy().sub(oldMousePos);
    fluid.addVelocity(
      Math.floor(mouse.relativePos.x * Fluid.N),
      Math.floor(mouse.relativePos.y * Fluid.N),
      diff.x,
      diff.y
    );
    oldMousePos = mouse.relativePos.copy();
  }
  ctx.fillRect(0, 0, canvas.width, canvas.height);

  fluid.step();
  fluid.renderD(ctx, canvas.width, canvas.height);
  fluid.fadeD();

  requestAnimationFrame(run);
}

paramConfig.onLoad(run);
