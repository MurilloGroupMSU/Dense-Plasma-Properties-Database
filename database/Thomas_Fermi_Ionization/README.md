
Atomic physics plays a crucial role in determining the properties of plasmas. A useful approximation is to assume that the electrons can be categorized as very strongly or very weakly interacting with the nuclei, which are often referred to as bound/free or core/valence, respectively. Knowing which electrons are in these cateogries is determined by temperature and pressure ionization, the determination of which is itelf a non-trivial calculation. There are three broad categories for performing such calculations:
1. Equilibrium, dilute plasmas: [Saha equation](https://en.wikipedia.org/wiki/Saha_ionization_equation) is used.
2. Non-equilibrium, dilute plasma: Rate equations, or their steady state solutions, are used; this is sometimes referred to as "atomic kinetics".
3. Dense, equilibrium plasmas: Pressure ionization is important, and average atom models are typically used. 

Here, an average-atom mean ionization state is computed using the [Thomas-Fermi electronic structure](https://en.wikipedia.org/wiki/Thomas%E2%80%93Fermi_model) model. The code is nearly trivial, owing to a fit given by Richard More.

Currently there are four versions of the code:
1. Single species Python (as a Jupyter notebook).
2. Single species Julia (as a Jupyter notebook).
3. Multi-species python (as a python module)
4. Multi-species C

If would like more information, or some really great papers to cite, consider reading/citing these:

[1] [Partial ionization in dense plasmas: Comparisons among average-atom density functional models](https://www.researchgate.net/publication/242376805_Partial_ionization_in_dense_plasmas_Comparisons_among_average-atom_density_functional_models), M. S. Murillo, Jon Weisheit, Stephanie B. Hansen, and M. W. C. Dharma-wardana, Phys. Rev. E 87, 063113 (2013).

[2] [Ionic transport in high-energy-density matter](https://www.researchgate.net/publication/300115955_Ionic_transport_in_high-energy-density_matter), Liam G. Stanton and Michael S. Murillo, Phys. Rev. E 93, (2016).
