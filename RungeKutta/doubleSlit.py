import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.animation as animation
from tqdm import tqdm

colors = [(1,1,1,c) for c in np.linspace(0,1,100)]
cmap_pot = mcolors.LinearSegmentedColormap.from_list("cmap_pot", colors, N = 5)

def RKstep(Hamiltonian, PsiT, length_x, length_y, k_now, k_next, scale_k):
    for j in range(length_y*length_x):
        for i in range(length_x*length_y):
            k_next[j] = k_next[j] + Hamiltonian[j * length_x*length_y + i] * (PsiT[i] + scale_k * k_now[i])
            

def setEdgesToZero(Psi, length_x, length_y):
    for j in range(length_y):
        for i in range(length_x):
            if i == 0 or i == length_x -1 or j == 0 or j == length_y -1: Psi[j * length_x + i] = 0

def getProbDens(Psi):
    newPsi = np.zeros(len(Psi))
    for i in range(len(Psi)):
        newPsi[i] = np.abs(Psi[i])
    return newPsi

#paramters
L = 10
dx = 0.1
dt = dx**2 / 4
Nx = int(L/dx) + 1
Ny = int(L/dx) + 1
N = Nx * Ny
Nt = 500
r = complex(0, 1/(dx**2))

x0 = L/5
y0 = L/2

kx = 15 * np.pi
sigma = 0.75

w = 0.2
s = 0.8
a = 0.4
V0 = 200

i0 = int(1 / (2 * dx) * (L - w))
i1 = int(1 / (2 * dx) * (L + w))

j0 = int(1 / (2 * dx) * (L + s) + a / dx)
j1 = int(1 / (2 * dx) * (L + s))
j2 = int(1 / (2 * dx) * (L - s))
j3 = int(1 / (2 * dx) * (L - s) - a / dx)

print("L=", L)
print("Nx=", Nx)
print("Ny=", Ny)
print("dx=", dx)
print("dt=", dt)

#initial wave packet
xarray = np.linspace(0, L, Nx)
yarray = np.linspace(0, L, Ny)
psi0 = np.zeros(Nx*Ny, dtype=complex)

for j in range(Ny):
    for i in range(Nx):
        if i == 0 or i == Nx -1 or j == 0 or j == Ny -1: psi0[j * Ny + i] = 0
        else:
            psi0[j * Ny + i] = np.exp(complex(0, kx * (xarray[i] - x0))) * \
            np.exp(-1.0 / (2 * sigma**2) * ((xarray[i] -x0)**2 + (yarray[j] - y0)**2))
print("wave packet initialized...")

x, y = np.meshgrid(xarray, yarray)
print("Meshgrid shape: ", x.shape)


"""
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(x, y, psi0.reshape(x.shape))
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
plt.show()
"""


#ptential
V = np.zeros(Nx*Ny, dtype=float)
for j in range(Ny):
    for i in range(Nx):
        if (i >= i0 and i <= i1) and ((j <= j3) or (j > j2 and j < j1) or (j >= j0)):
            V[j * Nx + i] = V0
print("potential set...")


#Hamiltonian
H = np.zeros(N*N, dtype=complex)
for k in range(N-1):
    j = math.floor(k/Nx)
    i = k % Nx

    #main diagonal
    H[k * N + k] = 4 * r - complex(0, V[j * Nx + i])
    #upper main diagonal
    if k+1 <= N: H[k * N + k + 1] = -r
    #lower main diagonal
    if k-1 >= 0: H[k * N + k - 1] = -r
    #upper lone diagonal
    if k+Nx <= N: H[k * N + k + Nx] = -r
    #lower lone diagonal
    if k-Nx >= 0: H[k * N + k - Nx]  = -r
print("Hamiltonian set...")



#run RungeKutta
psiT = psi0
psis = []
psi0 = getProbDens(psi0)
psis.append(psi0)
for t in tqdm(range(Nt)):
    """
    k0 = np.zeros(N, dtype=complex)
    k1 = np.zeros(N, dtype=complex)
    k2 = np.zeros(N, dtype=complex)
    k3 = np.zeros(N, dtype=complex)
    k4 = np.zeros(N, dtype=complex)

    RKstep(H, psiT, Nx, Ny, k0, k1, 0)
    RKstep(H, psiT, Nx, Ny, k1, k2, 0.5*dt)
    RKstep(H, psiT, Nx, Ny, k2, k3, 0.5*dt)
    RKstep(H, psiT, Nx, Ny, k3, k4, dt)
    """

    k1 = np.dot(H.reshape(N, N), psiT)
    k2 = np.dot(H.reshape(N, N), (psiT + 0.5 * dt * k1))
    k3 = np.dot(H.reshape(N, N), (psiT + 0.5 * dt * k2))
    k4 = np.dot(H.reshape(N, N), (psiT + 1.0 * dt * k3))


    for i in range(N):
        psiT[i] = psiT[i] + (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i])

    setEdgesToZero(psiT, Nx, Ny)
    newPsi = getProbDens(psiT)
    psis.append(newPsi)



psis = np.asarray(psis)

fig = plt.figure()
ax = fig.add_subplot(111, xlim=(0,L), ylim=(0,L))
ax.imshow(V.reshape(x.shape), extent=[0,L,0,L], cmap=cmap_pot, vmin=1, vmax=200, zorder=50, interpolation="none")
img = ax.imshow(psis[0].reshape(x.shape), extent=[0,L,0,L], cmap=plt.get_cmap("hot"), vmin=np.min(psis[0]), vmax=np.max(psis[0]), zorder=1, interpolation="none")
def animate(i):
    img.set_data(psis[i].reshape(x.shape))
    img.set_zorder(1)
anim=animation.FuncAnimation(fig, animate, interval=1, frames=np.arange(0, len(psis), 2), repeat=False, blit=0)
writer = animation.PillowWriter(fps=60)
anim.save("doubleSlit.gif", writer=writer)

#plt.show()

