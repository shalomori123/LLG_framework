# 2. Abstract

*(~250 words. Advisor feedback incorporated: no "weak effect" framing;
literature referenced generically.)*

**Parallel GPU Computing of Light–Magnetism Interaction**

When an electromagnetic wave propagates through a magnetic medium, the
magnetic component of the light — usually neglected in conventional optics —
couples to the magnetization of the material and produces a magneto-optical
response that can become significant in certain regimes, such as circularly
polarized light. Resolving this effect requires solving Maxwell's equations
and the Landau–Lifshitz–Gilbert (LLG) equation as a coupled, two-way system,
and at the relevant amplitudes double-precision arithmetic is needed to
preserve the signal. On a CPU such simulations take hours to days, which
motivates a parallel approach.

In this project we developed a GPU framework that integrates the coupled
Maxwell–LLG system with the FDTD method and a fourth-order Runge–Kutta LLG
step. A time step consists of two custom CUDA kernels, one for the magnetic
update (Faraday and LLG) and one for the electric update (Ampère), with one
GPU thread per grid cell. Neighbouring values are exchanged through the
cheapest level of the memory hierarchy that holds them: warp shuffles,
shared memory, then global memory. A single kernel source covers both single
and double precision, grids of any size, and batches of independent material
cases running in the same launches. Boundary conditions and the source run
inside the kernels without inter-block synchronization.

The implementation was validated against published results: reflection and
transmission agree to five decimal places, the magnetization magnitude is
conserved to 10⁻¹³, and the 10⁻⁹ magneto-optical signature is resolved. The
GPU framework runs 15×–180× faster than the CPU baseline, depending on grid
size.
