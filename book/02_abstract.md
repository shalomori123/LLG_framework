# 2. Abstract

*(~250 words. Advisor feedback incorporated: no "weak effect" framing;
literature referenced generically.)*

**Parallel GPU Computing of Light–Magnetism Interaction**

When an electromagnetic wave propagates through a magnetic medium, the
magnetic component of the light — usually neglected in conventional optics —
couples to the magnetization of the material and produces a magneto-optical
response that can become significant in certain regimes, such as circularly
polarized light. Resolving this effect requires solving Maxwell's equations
and the Landau–Lifshitz–Gilbert (LLG) equation as a coupled, two-way system
over millions of time steps. With the lab's previous CPU method such
simulations take hours to days, which motivates a parallel approach.

In this project we developed a GPU framework that integrates the coupled
Maxwell–LLG system with the finite-difference time-domain (FDTD) method and a
fourth-order Runge–Kutta LLG step. A time step consists of two custom CUDA
kernels, one for the magnetic update (Faraday's law and LLG) and one for the
electric update (Ampère's law), with one GPU thread per grid cell.
Neighbouring values are exchanged through the cheapest level of the memory
hierarchy that holds them: warp shuffles, shared memory, then global memory.
Boundary conditions and the source run inside the kernels without
inter-block synchronization, independent material cases share the same
launches, and the output is streamed to disk in bounded memory. A single
kernel source compiles to single or double precision, so the
accuracy–performance trade-off can be chosen per study.

The implementation reproduces the lab's previous reflection and transmission
results to five decimal places, conserves the magnetization magnitude to
10⁻¹³, and resolves the 10⁻⁹ magneto-optical signature. The GPU framework
runs 15×–180× faster than the CPU method, depending on grid size.
