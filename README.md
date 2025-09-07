Simulación Híbrida FEM-PDE de la Respuesta Sísmica en Edificios Altos: Un Enfoque Computacional en Python
Autor	Dr. Julio Rodrigo Baeza Pereyra, PhD
RESUMEN
El análisis de la respuesta sísmica de un edificio de varios pisos sometido a efectos de fuerzas sísmicas ha sido estudiado bajo diferentes enfoques de trabajo. En la actualidad se conocen bien dos métodos de estos: el de ecuaciones diferenciales parciales (PDE: Partial Diferencial Equations) y el de los elementos finitos (FEM: Finite Elements Method). En este trabajo se presenta un enfoque híbrido que combina los dos métodos en uno solo. Este enfoque se muestra como introducción a la aplicación de simulaciones de comportamiento estructural y es un inicio en los que tanto educadores como estudiantes de niveles de ingeniería pueden conocer las aplicaciones numéricas avanzadas en los computadores pc de escritorio o portátiles. Este trabajo solo se considera como uno de los que se puede partir hacia trabajos de tesis tanto de licenciaturas, como maestrías y si es posible, como antecedente a estudios de posgrado de niveles superiores en ingeniería estructural.
Palabras clave: Metodo de los Elementos Finitos, FEM, Ecuaciones Diferenciales Parciales, PDE, Métodos Acoplados, Enfoque Hybrido, Programación Python.
ABSTRACT
The analysis of the seismic response of a multi-storey building subjected to the effects of seismic forces has been studied under different work approaches. Two of these methods are currently well known: Partial Differential Equations (PDE) and Finite Elements Method (FEM). This paper presents a hybrid approach that combines the two methods into one. This approach is shown as an introduction to the application of structural behavior simulations and is a beginning in which both educators and students of engineering levels can learn about advanced numerical applications in desktop or laptop computers. This work is only considered as one of those that can be started towards thesis work for both bachelor's and master's degrees and, if possible, as a precursor to postgraduate studies at higher levels in structural engineering.
Keywords: Finite Element Method, FEM, Partial Differential Equations, PDE, Coupled Methods, Hybrid Approach, Python Programming.
MARCO TEÓRICO
Para la solución del problema de la respuesta estructural de edificios altos bajo los efectos de un sismo de duración t, se parte de la ecuación de equilibrio dinámico, la cual es una ecuación de derivadas parciales (PDE):
\rho\mathbit{A}\frac{\partial^\mathbf{2}\mathbit{u}\left(\mathbit{x},\mathbit{t}\right)}{\partial\mathbit{t}^\mathbf{2}}+\mathbit{C}\frac{\partial\mathbit{u}\left(\mathbit{x},\mathbit{t}\right)}{\partial\mathbit{t}}+\mathbit{EI}\frac{\partial^\mathbf{4}\mathbit{u}\left(\mathbit{x},\mathbit{t}\right)}{\partial\mathbit{x}^\mathbf{4}}=\mathbit{f}\left(\mathbit{x},\mathbit{t}\right)  	ecuación 1
Donde:
	es la densidad del material;
u(x,t): 	desplazamiento transversal en función de la altura x y el tiempo t
A: área de la sección transversal
c: coeficiente de amortiguamiento o de Rayleigh
E: módulo de elasticidad
I: momento de Inercia de la sección
f(x, t): fuerza externa (en este caso, sísmica aplicada a la base, que se transmite a los pisos superiores)

Discscretización (PDE a ODE)
\frac{\partial^\mathbf{4}\mathbit{u}}{\partial\mathbit{x}^\mathbf{4}}\approx\frac{\mathbit{u}_{\mathbit{i}-\mathbf{2}}-\mathbf{4}\mathbit{u}_{\mathbit{i}-\mathbf{1}}+\mathbf{6}\mathbit{u}_\mathbit{i}-\mathbf{4}\mathbit{u}_{\mathbit{i}+\mathbf{1}}+\mathbit{u}_{\mathbit{i}+\mathbf{2}}}{\mathbit{d}\mathbit{x}^\mathbf{4}} 	ecuación 2

Donde se discretiza esta PDE usando diferencias finitas para la derivada cuarta y esto genera la matriz de rigidez K_pde que representa el operador de flexión.


Acoplamiento FEM
El modelo FEM parte de la formulación clásica de dinámica estructural:
M\ddot{u}\left(t\right)+C\dot{u}\left(t\right)+Ku\left(t\right)=F\left(t\right) 	ecuación 3
Donde: 
M: matriz de masa
C: matriz de amortiguamiento
K: matriz de rigidez
u(t): vector de desplazamientos
v(t)=   vector de velocidades (la primera derivada de los desplazamientos)
 = vector de aceleraciones (la segunda derivada de los desplazamientos)
F(t): vector de fuerza

Luego se construye la matriz combinando masas, rigideces e inercias para un sistema discreto de elementos i y j, en este caso tipo pórtico, y luego las combina con la matriz del sistema de ecuaciones híbrido.
Sistema final acoplado. El sistema final que se puede resolver se define como:
M_{\mathrm{hybrid}}\ddot{u}\left(t\right)+C_{\mathrm{hybrid}}\dot{u}\left(t\right)+K_{\mathrm{hybrid}}u\left(t\right)=F\left(t\right)	ecuación 4

Después del acoplamiento, el sistema dinámico a resolver es:
\dot{u}\left(t\right)=v\left(t\right)	ecuación 5
\dot{v}\left(t\right)=M^{-1}\left(F\left(t\right)-Cv\left(t\right)-Ku\left(t\right)\right)	ecuación 6

Este sistema se puede representar de la siguiente manera:
\frac{d}{dt}\left[\begin{matrix}u\left(t\right)\\v\left(t\right)\\\end{matrix}\right]=\left[\begin{matrix}v\left(t\right)\\M^{-1}\left(F\left(t\right)-Cv\left(t\right)-Ku\left(t\right)\right)\\\end{matrix}\right]	ecuación 7
Lo que resulta en los siguientes vectores:
u\left(t\right)=\left[\begin{matrix}u_1\ \\u_2\\.\\.\\.\\u_n\\\end{matrix}\right] , v\left(t\right)=\ \left[\begin{matrix}v_1\\v_2\\.\\.\\.\\v_n\\\end{matrix}\right]	ecuación 8

Para cada piso i = 1, 2, …, n, la ecuación de velocidad es:
\dot{u_i}\left(t\right)=v_i\left(t\right) 	ecuación 9
\ddot{u_i}\left(t\right)=\frac{1}{m_i}\left(f_i\left(t\right)-\sum_{j=1}^{n}{c_{ij}\dot{u_j}\left(t\right)}-\sum_{j=1}^{n}{k_{ij}u_j\left(t\right)}\right) ecuación 10
Donde 
   = masa total del piso i (esto incluye la carga viva)
  = coeficiente de amortiguamiento entre pisos
  = coeficiente de rigidez entre pisos
  = fuerza sísmica en el piso i

