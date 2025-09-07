import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import diags
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

## ================= PARÁMETROS DEL MODELO HÍBRIDO ======
n_pisos = 100                   # Número de pisos
L_total = 3.5 * n_pisos         # Altura total (m)

E = 1.75e11                     # Módulo elástico del acero estructural (Pa)
rho = 6500.0                    # Densidad del acero (kg/m³)
A = 0.025                       # Área de sección transversal que componen los elementos, 
                                # en este caso constante para el edificio (m²)
I = 5.5e-4                      # Momento de inercia (m⁴)
alpha_pde = 0.03                # Amortiguamiento PDE
beta_fem = 0.003                # Amortiguamiento FEM

## ================= MODELO PDE =================
def setup_pde_system():
    dx = L_total / (n_pisos - 1)
    x = np.linspace(0, L_total, n_pisos)

    K_pde = (E * I / dx**4) * diags([1, -4, 6, -4, 1],
                                    [-2, -1, 0, 1, 2],
                                    shape=(n_pisos, n_pisos)).toarray()

    M_pde = rho * A * np.eye(n_pisos)
    C_pde = alpha_pde * M_pde + beta_fem * K_pde

    return K_pde, M_pde, C_pde

## ================= MODELO FEM =================
def setup_fem_system():
    n_nodos = n_pisos
    n_gdl = n_nodos * 3

    K_fem = np.zeros((n_gdl, n_gdl))
    M_fem = np.zeros((n_gdl, n_gdl))

    L_avg = L_total / (n_pisos - 1)

    def element_stiffness(E, A, I, L):
        k11 = E*A/L
        k22 = 12*E*I/L**3
        k33 = 4*E*I/L
        return np.array([
            [k11, 0, 0, -k11, 0, 0],
            [0, k22, 6*E*I/L**2, 0, -k22, 6*E*I/L**2],
            [0, 6*E*I/L**2, k33, 0, -6*E*I/L**2, 2*E*I/L],
            [-k11, 0, 0, k11, 0, 0],
            [0, -k22, -6*E*I/L**2, 0, k22, -6*E*I/L**2],
            [0, 6*E*I/L**2, 2*E*I/L, 0, -6*E*I/L**2, k33]
        ])
    for i in range(n_pisos - 1):
        idx = 3 * i
        K_elem = element_stiffness(E, A, I, L_avg)
        K_fem[idx:idx+6, idx:idx+6] += K_elem

    for i in range(n_pisos):
        M_fem[3*i, 3*i] = rho * A * L_avg / 2
        M_fem[3*i+1, 3*i+1] = rho * A * L_avg / 2

    C_fem = alpha_pde * M_fem + beta_fem * K_fem

    return K_fem, M_fem, C_fem

## ================= ACOPLAMIENTO HÍBRIDO =================
def coupled_system():
    K_pde, M_pde, C_pde = setup_pde_system()
    K_fem, M_fem, C_fem = setup_fem_system()

    hybrid_factor = 0.5

    K_hybrid = hybrid_factor * K_fem[::3, ::3] + (1-hybrid_factor) * K_pde
    M_hybrid = hybrid_factor * M_fem[::3, ::3] + (1-hybrid_factor) * M_pde
    C_hybrid = hybrid_factor * C_fem[::3, ::3] + (1-hybrid_factor) * C_pde

    return K_hybrid, M_hybrid, C_hybrid

## ================= EXCITACIÓN SÍSMICA =================
def sismo(t):
    return 1.5 * (np.sin(2*np.pi*0.25*t) +
                  0.4*np.sin(2*np.pi*1.2*t) +
                  0.2*np.sin(2*np.pi*3.0*t)) * np.exp(-0.05*t)

## ================= SISTEMA DINÁMICO =================
def dynamic_equations(t, y, K, M, C):
    n = len(K)
    u = y[:n]
    v = y[n:]

    F = -M[:, 0] * sismo(t)
    M_smooth = M + 1e-6 * np.eye(n)

    dudt = v
    dvdt = np.linalg.solve(M_smooth, F - C @ v - K @ u)

    return np.concatenate([dudt, dvdt])

## ================= SOLUCIÓN NUMÉRICA =================
K, M, C = coupled_system()
n = len(K)
y0 = np.zeros(2 * n)
t_eval = np.linspace(0, 30, 600)  # Más pasos para mayor resolución

sol = solve_ivp(dynamic_equations, [0, 30], y0, t_eval=t_eval,
                args=(K, M, C), method='BDF', rtol=1e-6)

u = sol.y[:n, :]

## ================= VISUALIZACIÓN =================
fig, ax = plt.subplots(figsize=(12, 8))
ax.set_xlim(-0.025, 0.025)
ax.set_ylim(0, L_total + 5)
ax.set_title("Simulación Híbrida PDE-FEM de Edificio de 100 Pisos", fontsize=14)
ax.set_xlabel("Desplazamiento (m)", fontsize=12)
ax.set_ylabel("Altura (m)", fontsize=12)
ax.grid(True)

height = np.linspace(0, L_total, n_pisos)
building, = ax.plot(u[:, 0], height, 'b-', lw=2, label='Edificio')
nodes = ax.scatter(u[::2, 0], height[::2], c='r', s=50, label='Nodos FEM')
ref_line = ax.axvline(x=0, color='k', linestyle='--', alpha=0.5)

def update(frame):
    scale = 10
    building.set_data(u[:, frame] * scale, height)
    nodes.set_offsets(np.column_stack((u[::2, frame] * scale, height[::2])))
    return building, nodes

ani = FuncAnimation(fig, update, frames=len(t_eval), blit=True, interval=50)
plt.legend()
plt.tight_layout()
plt.show()

