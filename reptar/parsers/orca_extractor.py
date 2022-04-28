# MIT License
# 
# Copyright (c) 2022, Alex M. Maldonado
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

from .extractor import extractor

class extractorORCA(extractor):
    """ORCA extractor

    Attributes
    ----------
    triggers : :obj:`tuple`
        A collection of triggers that activate the corresponding extractor.
        The trigger is a lambda function that returns True or False depending
        on the criteria and the name of the extractor method.
    parsed_info : :obj:`dict`
        Information parsed from files. Contains the following keys.

        ``system_info``
            Information specifying the system prior to any computation. Such
            as the initial cartesian coordinates, total system charge and
            multiplicity, etc.
        
        ``runtime_info``
            Contains information about setting up the job/calculation or running
            the job. Defining convergence criteria, parameters, etc.
        
        ``outputs``
            Results, requested or not, of the job. For example, SCF
            cycle values, optimized coordinates, trajectory, number of
            electrons, generated structures, etc.
    """
    def __init__(self):
        super().__init__()
    
    @property
    def triggers(self):
        trig = (
            (lambda line: True if ('Your calculation utilizes the auxiliary basis: ' in line.strip()) else False, 'aux_basis'),
            (lambda line: True if ('DFT GRID GENERATION' == line.strip() or 'Setting up the final grid:' == line.strip()) else False, 'grid_info'),
            (lambda line: True if ('CPCM SOLVATION MODEL' == line.strip()) else False, 'implicit_solvent'),
            (lambda line: True if ('TOTAL SCF ENERGY' == line.strip()) else False, 'scf_energies'),
            (lambda line: True if ('UHF SPIN CONTAMINATION' == line.strip()) else False, 'uhf_spin_contamination'),
            (lambda line: True if ('ORCA  MP2' == line.strip()) else False, 'mp_properties'),
            (lambda line: True if ('Wavefunction type' == line.strip()) else False, 'cc_method'),
            (lambda line: True if ('HF COUPLED CLUSTER ITERATIONS' == line.strip()[1:] or 'OPEN-SHELL COUPLED CLUSTER ITERATIONS' == line.strip()) else False, 'cc_properties'),
            (lambda line: True if ('SCF SETTINGS' == line.strip()) else False, 'scf_info'),
            (lambda line: True if ('MULLIKEN ATOMIC CHARGES' == line.strip()) else False, 'mulliken_charges'),
            (lambda line: True if ('LOEWDIN ATOMIC CHARGES' == line.strip()) else False, 'loewdin_charges'),
            (lambda line: True if ('DIPOLE MOMENT' == line.strip()) else False, 'dipole'),
            (lambda line: True if ('Geometry convergence' in line and '----------------------' in line) else False, 'geo_conv'),
            (lambda line: True if ('VIBRATIONAL FREQUENCIES' == line[:23]) else False, 'frequencies'),
            (lambda line: True if ('NORMAL MODES' == line[:12]) else False, 'normal_modes'),
            (lambda line: True if ('Temperature         ...' in line) else False, 'thermo_temp'),
            (lambda line: True if ('Zero point energy                ...' in line) else False, 'zpve'),
            (lambda line: True if ('Total thermal correction' in line) else False, 'thermal_corrections'),
            (lambda line: True if ('Thermal Enthalpy correction       ...' in line) else False, 'enthalpic_corr'),
            (lambda line: True if ('Final entropy term                ...' in line) else False, 'entropic_corr'),
        )
        return trig

    def aux_basis(self, f, line):
        """Information about auxiliary basis sets used in the calculation.

        Parameters
        ----------
        f : :obj:`io.TextIOWrapper`
            Buffered text stream of the file.
        line : :obj:`str`
            Parsed line from ``f``.
        
        Notes
        -----
        Example trigger text for this extractor.

        .. code-block:: text

            ----- AuxJ basis set information -----
            Your calculation utilizes the auxiliary basis: def2/J
        """
        line_split = line.strip().split(':')
        self.parsed_info['runtime_info']['basis_aux'] = line_split[-1].strip()

    def grid_info(self, f, line):
        """DFT integration grid information.

        All DFT calculations are performed with numerical integration, which
        means a integration grid must be specified. This is typically handeled
        with the Grid keyword. Furthermore, ORCA defaults to a multigrid
        approach where one grid is used for the SCF cycle, and a different
        (typically larger) grid is used for the final energy evaluation.

        After testing different keyword combinations, the Grid keywords are not
        always consistent with the final grid. Thus, we are going to directly
        parse the grid information from the output file instead of depending
        on the keywords.

        Parameters
        ----------
        f : :obj:`io.TextIOWrapper`
            Buffered text stream of the file.
        line : :obj:`str`
            Parsed line from ``f``.
        
        Notes
        -----
        Example trigger text for this extractor.

        .. code-block:: text

            -------------------
            DFT GRID GENERATION
            -------------------
        
        Here are the main parameters specifying the ORCA default grid levels.
        Default SCF grid is 2.

        +--------+---------------+----------+
        |  Grid  |  AngularGrid  |  IntAcc  |
        +========+===============+==========+
        |  1     |  Lebedev-50   |  4.34    |
        |  2     |  Lebedev-110  |  4.34    |
        |  3     |  Lebedev-194  |  4.34    |
        |  4     |  Lebedev-302  |  4.67    |
        |  5     |  Lebedev-434  |  5.01    |
        |  6     |  Lebedev-590  |  5.34    |
        |  7     |  Lebedev-770  |  5.67    |
        +--------+---------------+----------+
        """
        lebedev_to_level = {
            '50': 1, '110': 2, '194': 3, '302': 4, '434': 5, '590': 6, '770': 7
        }
        
        # SCF Cycle Grid
        if 'DFT GRID GENERATION' == line.strip():
            while 'Angular Grid (max. acc.)' not in line:
                line = next(f)
            lebedev_num = line.strip().split('-')[-1]
            self.parsed_info['runtime_info']['scf_grid_level'] = lebedev_to_level[lebedev_num]
        elif 'Setting up the final grid:' == line.strip():
            while 'Angular Grid (max. acc.)' not in line:
                line = next(f)
            lebedev_num = line.strip().split('-')[-1]
            self.parsed_info['runtime_info']['scf_grid_level_final'] = lebedev_to_level[lebedev_num]

    def implicit_solvent(self, f, line):
        """Implicit solvent properties.

        Parameters
        ----------
        f : :obj:`io.TextIOWrapper`
            Buffered text stream of the file.
        line : :obj:`str`
            Parsed line from ``f``.
        
        Notes
        -----
        Example trigger text for this extractor.

        .. code-block:: text

            --------------------
            CPCM SOLVATION MODEL
            --------------------
        """
        while 'Overall time for CPCM initialization' != line.strip()[:36]:
            # Solvent:                                     ... ACETONITRILE
            if 'Solvent:' == line[:8]:
                self.parsed_info['runtime_info']['solv_implicit_name'] = line.split()[2]
            
            line = next(f)

    
    def scf_energies(self, f, line):
        """The nuclear repulsion, one- and two-electron energy, and
        exchange-correlation energy after a SCF cycle.

        This is called directly after the ``'TOTAL SCF ENERGY'`` trigger, and 
        will terminate once the ``'SCF CONVERGENCE'`` trigger is reached.

        Instead of returning the energies themselves, we handle the creation and
        modification of ``scf_info`` here so any missing information (such as
        ``'energy_scf_xc'`` in MP2 calculations) is not an issue.

        Parameters
        ----------
        f : :obj:`io.TextIOWrapper`
            Buffered text stream of the file.
        line : :obj:`str`
            Parsed line from ``f``.
        
        Notes
        -----
        Example trigger text for this extractor.

        .. code-block:: text

            ----------------
            TOTAL SCF ENERGY
            ----------------
        """
        while 'SCF CONVERGENCE' != line.strip():
            # Nuclear Repulsion  :  135.87324654 Eh    3697.29901 eV
            if 'Nuclear Repulsion  :' in line:
                if 'energy_nuc_repul' not in self.parsed_info['outputs'].keys():
                    self.parsed_info['outputs']['energy_nuc_repul'] = []
                self.parsed_info['outputs']['energy_nuc_repul'].append(
                    float(line.split()[3])
                )

            # One Electron Energy: -674.26034691 Eh  -18347.55681 eV
            if 'One Electron Energy:' in line:
                if 'energy_scf_one_ele' not in self.parsed_info['outputs'].keys():
                    self.parsed_info['outputs']['energy_scf_one_ele'] = []
                self.parsed_info['outputs']['energy_scf_one_ele'].append(
                    float(line.split()[3])
                )

            # Two Electron Energy:  245.90403408 Eh    6691.38895 eV
            if 'Two Electron Energy:' in line:
                if 'energy_scf_two_ele' not in self.parsed_info['outputs'].keys():
                    self.parsed_info['outputs']['energy_scf_two_ele'] = []
                self.parsed_info['outputs']['energy_scf_two_ele'].append(
                    float(line.split()[3])
                )

            # E(XC)              :      -26.170411238000 Eh 
            if 'E(XC) ' in line:
                if 'energy_scf_xc' not in self.parsed_info['outputs'].keys():
                    self.parsed_info['outputs']['energy_scf_xc'] = []
                self.parsed_info['outputs']['energy_scf_xc'].append(
                    float(line.split()[2])
                )

            line = next(f)

    
    def uhf_spin_contamination(self, f, line):
        """Spin contamination from unrestricted Hartree-Fock calculations.

        Parameters
        ----------
        f : :obj:`io.TextIOWrapper`
            Buffered text stream of the file.
        line : :obj:`str`
            Parsed line from ``f``.
        
        Notes
        -----
        Example trigger text for this extractor.

        .. code-block:: text

            ----------------------
            UHF SPIN CONTAMINATION
            ----------------------
            
            Expectation value of <S**2>     :     0.750016
            Ideal value S*(S+1) for S=0.5   :     0.750000
            Deviation                       :     0.000016
        """

        for _ in range(0, 3):
            line = next(f)
        if line.strip() == 'Warning: in a DFT calculation there is little theoretical justification to':
            for _ in range(0, 4):
                line = next(f)
        
        if 'spin_sqrd_uhf_ideal' not in self.parsed_info['outputs'].keys():
            self.parsed_info['outputs']['spin_sqrd_uhf_ideal'] = []
        if 'spin_sqrd_uhf_calc' not in self.parsed_info['outputs'].keys():
            self.parsed_info['outputs']['spin_sqrd_uhf_calc'] = []

        # Expectation value of <S**2>     :     0.750016
        self.parsed_info['outputs']['spin_sqrd_uhf_calc'].append(
            float(line.split()[5])
        )
        line = next(f)
        # Ideal value S*(S+1) for S=0.5   :     0.750000
        self.parsed_info['outputs']['spin_sqrd_uhf_ideal'].append(
            float(line.split()[6])
        )

    
    def mp_properties(self, f, line):
        """Moller-Plesset calculation properties.

        This is called directly after the ``'ORCA  MP2 '`` trigger, and 
        will terminate once the ``'ORCA property calculations'`` trigger is reached.

        Instead of returning the energies themselves, we handle the creation and
        modification of ``mp_info`` here so any missing information is not an
        issue.

        Parameters
        ----------
        f : :obj:`io.TextIOWrapper`
            Buffered text stream of the file.
        line : :obj:`str`
            Parsed line from ``f``.
        
        Notes
        -----
        Example trigger text for this extractor.

        .. code-block:: text

            ----------------------------------------------------------
                                    ORCA  MP2 
            ----------------------------------------------------------
            
            Expectation value of <S**2>     :     0.750016
            Ideal value S*(S+1) for S=0.5   :     0.750000
            Deviation   
        """
        while '-     ORCA property calculations      *' != line.strip():
            # Freezing NCore=10 chemical core electrons
            if 'Freezing NCore=' == line.strip()[:15]:
                # if 'ele_frozen' not in self.data['keywords'].keys():
                #     self.data['outputs']['ele_frozen'] = []
                ele_frozen = int(line.split()[1][6:])
                self.parsed_info['outputs']['ele_frozen'] = ele_frozen

            #  MP2 CORRELATION ENERGY   :     -3.132364939 Eh
            if 'MP2 CORRELATION ENERGY' in line:
                if 'energy_correl_mp2' not in self.parsed_info['outputs'].keys():
                    self.parsed_info['outputs']['energy_correl_mp2'] = []
                if 'RI-MP2' in line:
                    index = 3
                else:
                    index = 4
                self.parsed_info['outputs']['energy_correl_mp2'].append(
                    float(line.split()[index])
                )
                break
            
            line = next(f)

    
    def cc_method(self, f, line):
        """Coupled cluster method information.

        Parameters
        ----------
        f : :obj:`io.TextIOWrapper`
            Buffered text stream of the file.
        line : :obj:`str`
            Parsed line from ``f``.
        
        Notes
        -----
        Example trigger text for this extractor.

        .. code-block:: text

            Wavefunction type
            -----------------
            Correlation treatment                      ...      CCSD     
            Single excitations                         ... ON 
        """
        while 'Number of correlated electrons' != line.strip()[:30]:
            # Frozen core treatment                      ... chemical core (0 el)
            # or
            # Frozen core treatment                      ... NO frozen core
            if 'Frozen core treatment' == line.strip()[:21]:
                if 'NO frozen core' == line.strip()[-14:]:
                    # self.parsed_info['runtime_info']['keywords']['frozen_core'] = False
                    pass
                else:
                    ele_frozen = int(line.split()[6][1:])
                    if ele_frozen == 0:
                        # self.parsed_info['runtime_info']['keywords']['frozen_core'] = False
                        pass
                    else:
                        # self.parsed_info['runtime_info']['keywords']['frozen_core'] = True
                        self.parsed_info['runtime_info']['ele_frozen'] = ele_frozen
            
            line = next(f)

    
    def cc_properties(self, f, line):
        """Coupled cluster properties.

        This is called directly after the ``'COUPLED CLUSTER ITERATIONS'``
        trigger, and will terminate once the ``'ORCA POPULATION ANALYSIS'``
        trigger is reached.

        Parameters
        ----------
        f : :obj:`io.TextIOWrapper`
            Buffered text stream of the file.
        line : :obj:`str`
            Parsed line from ``f``.
        
        Notes
        -----
        Example trigger text for this extractor.

        .. code-block:: text

            ------------------------------------------------
                              RHF COUPLED CLUSTER ITERATIONS
            ------------------------------------------------
        
        .. code-block:: text

            ------------------------------------------------
                              UHF COUPLED CLUSTER ITERATIONS
            ------------------------------------------------

        .. code-block:: text
        
            ----------------------------------------------------------
                                 OPEN-SHELL COUPLED CLUSTER ITERATIONS
            ----------------------------------------------------------
        """
        while 'ORCA POPULATION ANALYSIS' != line.strip() \
        and '*     ORCA property calculations      *' != line.strip():
            # Iter       E(tot)           E(Corr)          Delta-E          Residual     Time      <S|S>**1/2
            #  0     -7.469294707     -0.036553055      0.000000000      0.027013328    0.05      0.000000001
            #
            # or
            #
            # Iter       E(tot)           E(Corr)          Delta-E          Residual     Time
            #   0     -2.897443580     -0.035772194      0.000000000      0.027217829    0.00
            if line.strip()[:79] == 'Iter       E(tot)           E(Corr)          Delta-E          Residual     Time':
                # Extracts MP2 energies under the initial line.
                line = next(f)
                if 'energy_correl_mp2' not in self.parsed_info['outputs'].keys():
                    self.parsed_info['outputs']['energy_correl_mp2'] = []
                if 'mp2_total_energy' not in self.parsed_info['outputs'].keys():
                    self.parsed_info['outputs']['mp2_total_energy'] = []
                
                self.parsed_info['outputs']['energy_correl_mp2'].append(
                    float(line.split()[2])
                )
                self.parsed_info['outputs']['mp2_total_energy'].append(
                    float(line.split()[1])
                )
            
            # E(TOT)                                     ...     -7.473852176
            if line.strip()[:6] == 'E(TOT)':
                # Extracts total CCSD energy.
                if 'energy_ccsd' not in self.parsed_info['outputs'].keys():
                    self.parsed_info['outputs']['energy_ccsd'] = []
                
                self.parsed_info['outputs']['energy_ccsd'].append(
                    float(line.split()[2])
                )
            
            # T1 diagnostic                              ...      0.001316573 
            if line.strip()[:13] == 'T1 diagnostic':
                # Extracts T1 diagnostic.
                if 'diag_t1' not in self.parsed_info['outputs'].keys():
                    self.parsed_info['outputs']['diag_t1'] = []
                
                self.parsed_info['outputs']['diag_t1'].append(
                    float(line.split()[3])
                )
            
            # E(CCSD(T))                                 ...     -7.473882409
            if line.strip()[:10] == 'E(CCSD(T))':
                # Extracts total CCSD(T) energy..
                if 'energy_ccsd(t)' not in self.parsed_info['outputs'].keys():
                    self.parsed_info['outputs']['energy_ccsd(t)'] = []
                
                self.parsed_info['outputs']['energy_ccsd(t)'].append(
                    float(line.split()[2])
                )
            
            line = next(f)

    
    def scf_info(self, f, line):
        """Other scf information.

        This will be placed under the ``'keyword'`` JSON property.

        Parameters
        ----------
        f : :obj:`io.TextIOWrapper`
            Buffered text stream of the file.
        line : :obj:`str`
            Parsed line from ``f``.

        Returns
        -------
        :obj:`dict`
            Available SCF information that could contain the following keys.

            ``'approx_rij'``
                The resolution of identity (RI) approximation for the Coulomb
                (J) term.
            ``'approx_cosx'``
                The chain-of-spheres integration approximation to the exchange
                term (COSX).
        
        Notes
        -----
        Example trigger text for this extractor.

        .. code-block:: text

            ------------
            SCF SETTINGS
            ------------
        """
        while 'Total time needed     ' not in line.strip():
            # Hamiltonian:
            #  Ab initio Hamiltonian  Method          .... Hartree-Fock(GTOs)


            # General Settings:
            #  Integral files         IntName         .... al.chrg0.mult2-orca.sp.esp-ccsdt.anopvqz.vtightscf.sym-lambda0
            #  Hartree-Fock type      HFTyp           .... UHF
            if 'Ab initio Hamiltonian' == line.strip()[:21]:
                # We only include the HF type in the keywords.
                for _ in range(0, 5):
                    line = next(f)
                hf_type = line.split()[4]
                self.parsed_info['runtime_info']['hf_type'] = hf_type
            
            #  Number of Electrons    NEL             ....    3
            if 'Number of Electrons' == line.strip()[:19]:
                if 'n_ele' not in self.parsed_info['runtime_info'].keys():
                    self.parsed_info['runtime_info']['n_ele'] = []
                n_electrons = int(line.split()[5])
                self.parsed_info['runtime_info']['n_ele'].append(n_electrons)

            #  RI-approximation to the Coulomb term is turned on
            if 'RI-approximation to the Coulomb term is turned on' in line:
                self.parsed_info['runtime_info']['approx_rij'] = True

            #    RIJ-COSX (HFX calculated with COS-X)).... on
            if 'RIJ-COSX (HFX calculated with COS-X)' in line:
                self.parsed_info['runtime_info']['approx_cosx'] = True
            
            if 'RI-JK (J+K treated both via RI)' in line:
                self.parsed_info['runtime_info']['approx_rijk'] = True
            
            line = next(f)

    
    def mulliken_charges(self, f, line):
        """Mulliken atomic charges in same order as atomic coordinates.

        Parameters
        ----------
        f : :obj:`io.TextIOWrapper`
            Buffered text stream of the file.
        line : :obj:`str`
            Parsed line from ``f``.
        
        Notes
        -----
        Example trigger text for this extractor.

        .. code-block:: text

            -----------------------
            MULLIKEN ATOMIC CHARGES
            -----------------------
        """
        line = self.skip_lines(f, 2)

        # Creates initial mulliken_charges property.
        if 'mulliken_charges' not in self.parsed_info['outputs'].keys():
            # self.parsed_info['outputs']['mulliken_charges'] = []
            pass
        
        # Appends Mulliken charges to a new item for every structure.
        # self.parsed_info['outputs']['mulliken_charges'].append([])
        while 'Sum of atomic charges' not in line:
            line_split = line.split(':')
            """
            self.parsed_info['outputs']['mulliken_charges'][-1].append(
                float(line_split[-1])
            )
            """
            line = next(f)

    
    def loewdin_charges(self, f, line):
        """Loewdin atomic charges in same order as atomic coordinates.

        Parameters
        ----------
        f : :obj:`io.TextIOWrapper`
            Buffered text stream of the file.
        line : :obj:`str`
            Parsed line from ``f``.
        
        Notes
        -----
        Example trigger text for this extractor.

        .. code-block:: text

            ----------------------
            LOEWDIN ATOMIC CHARGES
            ----------------------
        """
        line = self.skip_lines(f, 2)

        # Creates initial loewdin_charges property.
        if 'loewdin_charges' not in self.parsed_info['outputs'].keys():
            # self.parsed_info['outputs']['loewdin_charges'] = []
            pass
        
        # Appends Loewdin charges to a new item for every structure.
        # self.parsed_info['outputs']['loewdin_charges'].append([])
        while '' != line.strip():
            line_split = line.split(':')
            """
            self.parsed_info['outputs']['loewdin_charges'][-1].append(
                float(line_split[-1])
            )
            """
            line = next(f)

    
    def dipole(self, f, line):
        """The X, Y, and Z dipole components.

        Final QCJSON specifies the method of the dipole moment (e.g.,
        ``'scf_dipole_moment'``, ``'mp2_dipole_moment'``). For now, we just
        store it as ``'dipole_moment'``.

        Parameters
        ----------
        f : :obj:`io.TextIOWrapper`
            Buffered text stream of the file.
        line : :obj:`str`
            Parsed line from ``f``.
        
        Notes
        -----
        Example trigger text for this extractor.

        .. code-block:: text

            -------------
            DIPOLE MOMENT
            -------------
        """
        if 'dipole_moment' not in self.parsed_info['outputs'].keys():
            # self.parsed_info['outputs']['dipole_moment'] = []
            pass
        
        while 'Total Dipole Moment    :' not in line:
            line = next(f)
        line_split = line.split()
        dipole = [float(line_split[4]), float(line_split[5]), float(line_split[6])]
        #self.parsed_info['outputs']['dipole_moment'].append(dipole)

    
    def _add_geo_conv(self, info_label, line):
        """Parse and add geometric convergence info to data.

        Parameters
        ----------
        info_label : :obj:`str`
            Label for geometric convergence criteria.
        line : :obj:`str`
            Line from output file to extract information from.
        """
        split_line = line.split()
        value = float(split_line[2])
        target = float(split_line[3])
        if f'geo_{info_label}_target' not in self.parsed_info['runtime_info'].keys():
            self.parsed_info['runtime_info'][f'geo_{info_label}_target'] = target
        try:
            self.parsed_info['runtime_info'][f'geo_{info_label}_value'].append(value)
        except KeyError:
            self.parsed_info['runtime_info'][f'geo_{info_label}_value'] = [value]

    
    def geo_conv(self, f, line):
        """Extract geometric convergence values and tolerance.

        Parameters
        ----------
        f : :obj:`io.TextIOWrapper`
            Buffered text stream of the file.
        line : :obj:`str`
            Parsed line from ``f``.
        
        Notes
        -----
        Example trigger text for this extractor.

        .. code-block:: text

                                  .--------------------.
            ----------------------|Geometry convergence|-------------------------
            Item                value                   Tolerance       Converged
            ---------------------------------------------------------------------
        """
        while 'Max(Dihed)' not in line and 'Max(Improp)' not in line:
            if 'Energy change' in line:
                self._add_geo_conv('energy_change', line)
            elif 'RMS gradient' in line:
                self._add_geo_conv('rms_gradient', line)
            elif 'MAX gradient' in line:
                self._add_geo_conv('max_gradient', line)
            elif 'RMS step' in line:
                self._add_geo_conv('rms_step', line)
            elif 'MAX step' in line:
                self._add_geo_conv('max_step', line)
            line = next(f)

    
    def frequencies(self, f, line):
        """Extract vibrational frequencies, freq_vibs. Includes 0.00 frequencies.

        Based on https://github.com/MolSSI/QCSchema/pull/50#issuecomment-499155251.

        Parameters
        ----------
        f : :obj:`io.TextIOWrapper`
            Buffered text stream of the file.
        line : :obj:`str`
            Parsed line from ``f``.
        
        Notes
        -----
        Example trigger text for this extractor.

        .. code-block:: text

            -----------------------
            VIBRATIONAL FREQUENCIES
            -----------------------
        """
        # Skips the following lines.
        # -----------------------
        # 
        # Scaling factor for frequencies =  1.000000000
        # 
        #    0:         0.00 cm**-1
        line = self.skip_lines(f, 5)

        # Sets up data.
        if 'freq_vib' not in self.parsed_info['outputs'].keys():
            self.parsed_info['outputs']['freq_vib'] = []

        vibfreqs = []
        while line.strip() != '':
            freq = float(line.split()[1])
            vibfreqs.append(freq)
            
            line = next(f)

        self.parsed_info['outputs']['freq_vib'].append(vibfreqs)

    
    def normal_modes(self, f, line):
        """Extract normalized, mass-weighted normal modes, q.

        Based on https://github.com/MolSSI/QCSchema/pull/50#issuecomment-499155251.

        Parameters
        ----------
        f : :obj:`io.TextIOWrapper`
            Buffered text stream of the file.
        line : :obj:`str`
            Parsed line from ``f``.
        
        Notes
        -----
        Example trigger text for this extractor.

        .. code-block:: text

            ------------
            NORMAL MODES
            ------------
        """
        line = self.skip_lines(f, 8)

        # Sets up data.
        if 'q' not in self.parsed_info['outputs'].keys():
            self.parsed_info['outputs']['normal_modes'] = []

        q = []
        mode = 0
        while line.strip() != '-----------':
            if '.' not in line:
                # Skips lines that print column numbers. For example,
                #   0          1          2          3          4          5 
                # Also restarts mode index.
                mode = 0
            else:
                disp = [float(i) for i in line.split()[1:]]
                try:
                    q[mode].extend(disp)
                except IndexError:
                    q.append(disp)
                mode += 1
            line = next(f)
        self.parsed_info['outputs']['normal_modes'].append(q)

    
    def thermo_temp(self, f, line):
        """Extracts temperature used for thermochemistry.

        Parameters
        ----------
        f : :obj:`io.TextIOWrapper`
            Buffered text stream of the file.
        line : :obj:`str`
            Parsed line from ``f``.
        
        Notes
        -----
        Example trigger text for this extractor.

        .. code-block:: text

            Temperature         ... 298.15 K
        """
        # Sets up data.
        if 'temp_thermochem' not in self.parsed_info['outputs'].keys():
            self.parsed_info['outputs']['temp_thermochem'] = []
        temp = float(line.split()[2])
        self.parsed_info['outputs']['temp_thermochem'].append(temp)
        # Needs to move off this line so extract can continue.
        line = next(f)

    
    def zpve(self, f, line):
        """Extract zero-point vibrational energy correction.

        Parameters
        ----------
        f : :obj:`io.TextIOWrapper`
            Buffered text stream of the file.
        line : :obj:`str`
            Parsed line from ``f``.
        
        Notes
        -----
        Example trigger text for this extractor.

        .. code-block:: text

            Zero point energy                ...      0.06455884 Eh      40.51 kcal/mol
        """
        if 'correc_zpe' not in self.parsed_info['outputs'].keys():
            self.parsed_info['outputs']['correc_zpe'] = []
        zpve = float(line.split()[4])
        self.parsed_info['outputs']['correc_zpe'].append(zpve)
        # Needs to move off this line so extract can continue.
        line = next(f)

    
    def thermal_corrections(self, f, line):
        """Extracts thermal vibrational, rotational, and translational
        corrections.

        Parameters
        ----------
        f : :obj:`io.TextIOWrapper`
            Buffered text stream of the file.
        line : :obj:`str`
            Parsed line from ``f``.
        
        Notes
        -----
        Example trigger text for this extractor.

        .. code-block:: text

            Total thermal correction                  0.00868494 Eh       5.45 kcal/mol
        """
        if 'correc_thermal' not in self.parsed_info['outputs'].keys():
            self.parsed_info['outputs']['correc_thermal'] = []
        thermal = float(line.split()[3])
        self.parsed_info['outputs']['correc_thermal'].append(thermal)
        line = next(f)  # Needs to move off this line so extract can continue.

    
    def enthalpic_corr(self, f, line):
        """Extract enthalpic corrections; the difference between H_298 and
        E_298.

        Parameters
        ----------
        f : :obj:`io.TextIOWrapper`
            Buffered text stream of the file.
        line : :obj:`str`
            Parsed line from ``f``.
        
        Notes
        -----
        Example trigger text for this extractor.

        .. code-block:: text

            Thermal Enthalpy correction       ...      0.00094421 Eh       0.59 kcal/mol
        """
        if 'correc_enthalpy' not in self.parsed_info['outputs'].keys():
            self.parsed_info['outputs']['correc_enthalpy'] = []
        enthalpy = float(line.split()[4])
        self.parsed_info['outputs']['correc_enthalpy'].append(enthalpy)
        line = next(f)  # Needs to move off this line so extract can continue.

    
    def entropic_corr(self, f, line):
        """Extract translational, rotational, and vibrational entropic
        corrections to enthalpy for Gibbs free energy (i.e., T*S_298).

        Parameters
        ----------
        f : :obj:`io.TextIOWrapper`
            Buffered text stream of the file.
        line : :obj:`str`
            Parsed line from ``f``.
        
        Notes
        -----
        Example trigger text for this extractor.
        
        .. code-block:: text

            Final entropy term                ...      0.04136711
        """
        if 'correc_entropy' not in self.parsed_info['outputs'].keys():
            self.parsed_info['outputs']['correc_entropy'] = []
        entropy = float(line.split()[4])
        self.parsed_info['outputs']['correc_entropy'].append(entropy)
        # Needs to move off this line so extract can continue.
        line = next(f)

