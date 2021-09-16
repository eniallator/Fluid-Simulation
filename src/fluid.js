class Fluid {
  static #N = 64;
  static #iter = 4;

  #dt;
  #diff;
  #visc;

  #s;
  #density;
  #Vx;
  #Vy;
  #Vx0;
  #Vy0;

  static get N() {
    return this.#N;
  }

  constructor(dt, diffusion, viscosity) {
    this.#dt = dt;
    this.#diff = diffusion;
    this.#visc = viscosity;

    this.#s = new Array(Fluid.#N * Fluid.#N).fill(0.0);
    this.#density = new Array(Fluid.#N * Fluid.#N).fill(0.0);
    this.#Vx = new Array(Fluid.#N * Fluid.#N).fill(0.0);
    this.#Vy = new Array(Fluid.#N * Fluid.#N).fill(0.0);
    this.#Vx0 = new Array(Fluid.#N * Fluid.#N).fill(0.0);
    this.#Vy0 = new Array(Fluid.#N * Fluid.#N).fill(0.0);
  }

  #IX(x, y) {
    return clamp(x, 0, Fluid.#N - 1) + clamp(y, 0, Fluid.#N - 1) * Fluid.#N;
  }

  addDensity(x, y, amount) {
    const index = this.#IX(x, y);

    this.#density[index] += amount;
  }

  addVelocity(x, y, amountX, amountY) {
    const index = this.#IX(x, y);

    this.#Vx[index] += amountX;
    this.#Vy[index] += amountY;
  }

  diffuse(b, x, x0, diff, dt) {
    const a = dt * diff * (Fluid.#N - 2) * (Fluid.#N - 2);
    this.#lin_solve(b, x, x0, a, 1 + 4 * a);
  }

  #lin_solve(b, x, x0, a, c) {
    const cRecip = 1 / c;
    for (let k = 0; k < Fluid.#iter; k++) {
      for (let j = 1; j < Fluid.#N - 1; j++) {
        for (let i = 1; i < Fluid.#N - 1; i++) {
          x[this.#IX(i, j)] =
            (x0[this.#IX(i, j)] +
              a *
                (x[this.#IX(i + 1, j)] +
                  x[this.#IX(i - 1, j)] +
                  x[this.#IX(i, j + 1)] +
                  x[this.#IX(i, j - 1)])) *
            cRecip;
        }
      }

      this.#set_bnd(b, x);
    }
  }

  project(velocX, velocY, p, div) {
    for (let j = 1; j < Fluid.#N - 1; j++) {
      for (let i = 1; i < Fluid.#N - 1; i++) {
        div[this.#IX(i, j)] =
          (-0.5 *
            (velocX[this.#IX(i + 1, j)] -
              velocX[this.#IX(i - 1, j)] +
              velocY[this.#IX(i, j + 1)] -
              velocY[this.#IX(i, j - 1)])) /
          Fluid.#N;
        p[this.#IX(i, j)] = 0;
      }
    }

    this.#set_bnd(0, div);
    this.#set_bnd(0, p);
    this.#lin_solve(0, p, div, 1, 4);

    for (let j = 1; j < Fluid.#N - 1; j++) {
      for (let i = 1; i < Fluid.#N - 1; i++) {
        velocX[this.#IX(i, j)] -=
          0.5 * (p[this.#IX(i + 1, j)] - p[this.#IX(i - 1, j)]) * Fluid.#N;
        velocY[this.#IX(i, j)] -=
          0.5 * (p[this.#IX(i, j + 1)] - p[this.#IX(i, j - 1)]) * Fluid.#N;
      }
    }

    this.#set_bnd(1, velocX);
    this.#set_bnd(2, velocY);
  }

  advect(b, d, d0, velocX, velocY, dt) {
    let i0, i1, j0, j1;

    let dtx = dt * (Fluid.#N - 2);
    let dty = dt * (Fluid.#N - 2);

    let s0, s1, t0, t1;
    let x, y;

    for (let j = 1; j < Fluid.#N - 1; j++) {
      for (let i = 1; i < Fluid.#N - 1; i++) {
        x = clamp(i - dtx * velocX[this.#IX(i, j)], 0.5, Fluid.#N + 0.5);
        y = clamp(j - dty * velocY[this.#IX(i, j)], 0.5, Fluid.#N + 0.5);

        i0 = Math.floor(x);
        i1 = i0 + 1;
        j0 = Math.floor(y);
        j1 = j0 + 1;

        s1 = x - i0;
        s0 = 1 - s1;
        t1 = y - j0;
        t0 = 1 - t1;

        d[this.#IX(i, j)] =
          s0 * (t0 * d0[this.#IX(i0, j0)] + t1 * d0[this.#IX(i0, j1)]) +
          s1 * (t0 * d0[this.#IX(i1, j0)] + t1 * d0[this.#IX(i1, j1)]);
      }
    }

    this.#set_bnd(b, d);
  }

  #set_bnd(b, x) {
    const N = Fluid.#N;
    for (let i = 1; i < N - 1; i++) {
      x[this.#IX(i, 0)] = b == 2 ? -x[this.#IX(i, 1)] : x[this.#IX(i, 1)];
      x[this.#IX(i, N - 1)] =
        b == 2 ? -x[this.#IX(i, N - 2)] : x[this.#IX(i, N - 2)];
    }
    for (let j = 1; j < N - 1; j++) {
      x[this.#IX(0, j)] = b == 1 ? -x[this.#IX(1, j)] : x[this.#IX(1, j)];
      x[this.#IX(N - 1, j)] =
        b == 1 ? -x[this.#IX(N - 2, j)] : x[this.#IX(N - 2, j)];
    }

    x[this.#IX(0, 0)] = 0.5 * (x[this.#IX(1, 0)] + x[this.#IX(0, 1)]);
    x[this.#IX(0, N - 1)] =
      0.5 * (x[this.#IX(1, N - 1)] + x[this.#IX(0, N - 2)]);
    x[this.#IX(N - 1, 0)] =
      0.5 * (x[this.#IX(N - 2, 0)] + x[this.#IX(N - 1, 1)]);
    x[this.#IX(N - 1, N - 1)] =
      0.5 * (x[this.#IX(N - 2, N - 1)] + x[this.#IX(N - 1, N - 2)]);
  }

  step() {
    const visc = this.#visc;
    const diff = this.#diff;
    const dt = this.#dt;
    const Vx = this.#Vx;
    const Vy = this.#Vy;
    const Vx0 = this.#Vx0;
    const Vy0 = this.#Vy0;
    const s = this.#s;
    const density = this.#density;

    this.diffuse(1, Vx0, Vx, visc, dt);
    this.diffuse(2, Vy0, Vy, visc, dt);

    this.project(Vx0, Vy0, Vx, Vy);

    this.advect(1, Vx, Vx0, Vx0, Vy0, dt);
    this.advect(2, Vy, Vy0, Vx0, Vy0, dt);

    this.project(Vx, Vy, Vx0, Vy0);

    this.diffuse(0, s, density, diff, dt);
    this.advect(0, density, s, Vx, Vy, dt);
  }

  renderD(ctx, width, height) {
    const cellWidth = Math.ceil(width / Fluid.#N);
    const cellHeight = Math.ceil(height / Fluid.#N);
    for (let i = 0; i < Fluid.#N; i++) {
      for (let j = 0; j < Fluid.#N; j++) {
        const shade = 255 * this.#density[this.#IX(i, j)];
        ctx.fillStyle = `rgb(${shade}, ${shade}, ${shade})`;
        ctx.fillRect(i * cellWidth, j * cellHeight, cellWidth, cellHeight);
      }
    }
  }

  fadeD() {
    for (let i = 0; i < this.#density.length; i++) {
      this.#density[i] = clamp(this.#density[i] - 1, 0, 255);
    }
  }
}

TimeAnalysis.registerMethods(Fluid);
