# MIT License
#
# Copyright (c) 2022-2023, Alex M. Maldonado
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

# pylint: disable=line-too-long

from .extractor import Extractor
from ..logger import ReptarLogger

log = ReptarLogger(__name__)


class ExtractorXTB(Extractor):
    r"""xTB extractor"""

    # We always pass the file object into methods here.
    # pylint: disable=unused-argument

    @property
    def triggers(self):
        trig = (
            (lambda line: bool("* xtb version" in line), "xtb_version"),
            (
                lambda line: bool("program call               :" in line.strip()),
                "run_type",
            ),
            (
                lambda line: bool("number of electrons        :" in line.strip()),
                "n_electrons",
            ),
            (
                lambda line: bool("charge                     :" in line.strip()),
                "charge",
            ),
            (
                lambda line: bool("spin                       :" in line.strip()),
                "multiplicity",
            ),
            (lambda line: bool("> wall" == line.strip()), "wall_pot"),
            (
                lambda line: (
                    ":                      SETUP                      :"
                    == line.strip()
                ),
                "gfn_setup",
            ),
            (
                lambda line: bool(
                    "::                     SUMMARY                     ::"
                    in line.strip()
                ),
                "summary_energies",
            ),
            (
                lambda line: bool(
                    "|               Molecular Dynamics                |"
                    in line.strip()
                ),
                "md_setup",
            ),
            (
                lambda line: bool(
                    "A N C O P T" == line.strip("| \n")
                    or "L-ANC optimizer" == line.strip("| \n")
                ),
                "opt_data",
            ),
            (
                lambda line: bool("average properties" == line.strip()),
                "md_avg_props",
            ),
            (
                lambda line: bool("thermostating problem" == line.strip()),
                "thermostat_prob",
            ),
            (
                lambda line: bool("normal termination of xtb" == line.strip()),
                "success",
            ),
        )
        return trig

    def xtb_version(self, f, line):
        r"""Version of xtb.

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

            * xtb version 6.3.3 (5b13467) compiled by 'ehlert@majestix' on 2020-09-17
        """
        line_split = line.split()
        version = line_split[3]
        self.parsed_info["runtime_info"]["prov_version"] = version
        next(f)

    def run_type(self, f, line):
        r"""Calculation type (e.g., opt, sp, grad).

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

            program call               : xtb 5h2o.example.xyz --scc --gfn 2 --charge 0
        """
        if "--scc" in line:
            driver = "energy"
        elif "--opt" in line:
            driver = "optimization"
        elif "--grad" in line:
            driver = "gradient"
        elif "--ohess" in line or "--hess" in line:
            driver = "frequency"
        elif "--omd" in line or "--md" in line:
            driver = "molecular_dynamics"
        self.parsed_info["runtime_info"]["calc_driver"] = driver
        next(f)

    def n_electrons(self, f, line):
        r"""Total number of electrons.

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

            number of electrons        :                  1688
        """
        n_ele = line.strip(" :").split()[-1]
        self.parsed_info["system_info"]["n_ele"] = int(n_ele)

    def charge(self, f, line):
        r"""Overall system charge.

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

            charge                     :                     0
        """
        _, _, charge = line.strip(" :").split()
        self.parsed_info["system_info"]["charge"] = int(charge)

    def multiplicity(self, f, line):
        r"""Overall system multiplicity.

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

            spin                       :                   0.0
        """
        _, _, spin = line.strip().split()
        mult = 2 * float(spin) + 1
        self.parsed_info["system_info"]["mult"] = int(mult)

    def wall_pot(self, f, line):
        r"""Wall potentials.

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

            > wall
            -> potential = logfermi
            -> sphere: 23.62158, all
            --> 1: 23.62158
            --> 2: all
        """
        if "wall_potential" in self.parsed_info:
            pots = self.parsed_info["runtime_info"]["wall_potential"]
        else:
            pots = []
        idx_start = len(pots)

        line = next(f)
        pot_type = line.split()[3]

        line = next(f)
        while "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" not in line:
            if "-> sphere:" == line[:10]:
                shape_type = line.split(":")[0].split(" ")[1]
                line = self.skip_lines(f, 2)
                atoms_constrained = ""
                while "spherical wallpotenial" != line[:22]:
                    line_split = line.strip().split()
                    if "all" == line_split[-1]:
                        atoms_constrained = "all"
                    else:
                        atoms_constrained += f"{line_split[-1]},"
                    line = next(f)
                if atoms_constrained[-1] == ",":
                    atoms_constrained = atoms_constrained[:-1]
                sphere_radius = float(line.split()[-2])
                pots.append(
                    {
                        "pot_type": pot_type,
                        "shape": shape_type,
                        "sphere_radius": sphere_radius,
                        "atoms_constrained": atoms_constrained,
                    }
                )
            elif "-> ellipsoid:" == line[:13]:
                pass

            if "-> temp =" == line[:9]:
                temp = float(line.split()[-1])
                line = next(f)
                beta = float(line.split()[-1])
                for i in range(idx_start, len(pots)):
                    pots[i]["logfermi_temp"] = temp
                    pots[i]["logfermi_beta"] = beta

            line = next(f)

        self.parsed_info["runtime_info"]["wall_potential"] = pots

    def gfn_setup(self, f, line):
        r"""GFN calculation setup.

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

            ...................................................
            :                      SETUP                      :
            :.................................................:
        """
        line = self.skip_lines(f, 2)
        while ".............................." not in line.strip():

            # :  # basis functions                1446          :
            if " # basis functions" in line:
                self.parsed_info["runtime_info"]["basis_n_func"] = int(
                    line.strip(" :").split()[-2]
                )

            # :  Hamiltonian                  GFN2-xTB          :
            if ":  Hamiltonian " in line:
                self.parsed_info["runtime_info"]["hamiltonian"] = line.strip(
                    " :"
                ).split()[1]

            line = next(f)

    def summary_energies(self, f, line):
        r"""Extracts energies listed in SUMMARY box.

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

            :::::::::::::::::::::::::::::::::::::::::::::::::::::
            ::                     SUMMARY                     ::
            :::::::::::::::::::::::::::::::::::::::::::::::::::::
            :: total energy             -25.404338997931 Eh    ::
            :: gradient norm              0.048763527905 Eh/a0 ::
            :: HOMO-LUMO gap             12.772900340050 eV    ::
            ::.................................................::
            :: SCC energy               -25.580127846040 Eh    ::
            :: -> isotropic ES            0.153316435656 Eh    ::
            :: -> anisotropic ES         -0.020934020006 Eh    ::
            :: -> anisotropic XC         -0.003785806991 Eh    ::
            :: -> dispersion             -0.004190013801 Eh    ::
            :: repulsion energy           0.175746000529 Eh    ::
            :: add. restraining           0.000000000000 Eh    ::
            :: total charge              -0.000000000000 e     ::
            :::::::::::::::::::::::::::::::::::::::::::::::::::::
        """
        if self.parsed_info["runtime_info"]["calc_driver"] != "molecular_dynamics":
            if "energy_scf" not in self.parsed_info["outputs"]:
                self.parsed_info["outputs"]["energy_scf"] = []

            if "energy_nuc_repul" not in self.parsed_info["outputs"]:
                self.parsed_info["outputs"]["energy_nuc_repul"] = []

            while "" != line.strip():
                if "total energy" in line:
                    line_split = line.split()
                    energy_scf = float(line_split[3])

                if "repulsion energy" in line:
                    line_split = line.split()
                    energy_nuc_repul = float(line_split[3])

                line = next(f)

            self.parsed_info["outputs"]["energy_scf"].append(energy_scf)
            self.parsed_info["outputs"]["energy_nuc_repul"].append(energy_nuc_repul)
        # For MD simulations.
        else:
            if "energy_pot" not in self.parsed_info["outputs"]:
                self.parsed_info["outputs"]["energy_pot"] = []

            while "" != line.strip():
                if "total energy" in line:
                    line_split = line.split()
                    energy_pot = float(line_split[3])
                    self.parsed_info["outputs"]["energy_pot"].append(energy_pot)

                line = next(f)

    def md_setup(self, f, line):
        r"""Extracts MD

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

                     -------------------------------------------------
                    |               Molecular Dynamics                |
                     -------------------------------------------------
            trajectories on xtb.trj or xtb.trj.<n>

            MD time /ps        :    5.00
            dt /fs             :    1.00
            SCC accuracy       :    1.00
            temperature /K     :  500.00
            max steps          :  5000
        """
        while "MD time /ps" not in line:
            line = next(f)

        #  MD time /ps        :    5.00
        self.parsed_info["runtime_info"]["md_duration"] = float(line.split()[4])

        line = next(f)
        #  dt /fs             :    1.00
        self.parsed_info["runtime_info"]["t_step"] = float(line.split()[3])

        line = next(f)
        #  SCC accuracy       :    1.00
        self.parsed_info["runtime_info"]["xtb_scc_accuracy"] = float(line.split()[3])

        line = next(f)
        #  temperature /K     :  500.00
        self.parsed_info["runtime_info"]["thermostat_temp"] = float(line.split()[3])

        while "dumpstep(trj) /fs" not in line:
            line = next(f)

        #  dumpstep(trj) /fs  :    5.00     5
        self.parsed_info["runtime_info"]["md_steps_dump_traj"] = float(line.split()[4])

        while "time (ps)    <Epot>" not in line:

            #  H atoms mass (amu) :     1
            if "H atoms mass (amu) :" == line.strip()[:20]:
                self.parsed_info["runtime_info"]["mass_h"] = float(line.split()[-1])

            # TODO: Parse shake on information.
            if "SHAKE off" == line.strip():
                self.parsed_info["runtime_info"]["algo_shake"] = False

            if "Berendsen THERMOSTAT on" == line.strip():
                self.parsed_info["runtime_info"]["thermostat_type"] = "Berendsen"

            if "RESTART" == line.strip():
                self.parsed_info["runtime_info"]["md_restarted"] = True

            line = next(f)

        if "md_restarted" not in self.parsed_info["runtime_info"]:
            self.parsed_info["runtime_info"]["md_restarted"] = False

    def opt_data(self, f, line):
        r"""All information during optimization routine. Incldues setup, energy.

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

             -----------------------------------------------------------
            |                   =====================                   |
            |                        A N C O P T                        |
            |                   =====================                   |
            |               Approximate Normal Coordinate               |
            |                Rational Function Optimizer                |
             -----------------------------------------------------------

        .. code-block:: text

             -----------------------------------------------------------
            |                       L-ANC optimizer                     |
             -----------------------------------------------------------
        """

        # The following lines signal the termination of the optimization routine.
        # *** GEOMETRY OPTIMIZATION CONVERGED AFTER 73 ITERATIONS ***
        # *** FAILED TO CONVERGE GEOMETRY OPTIMIZATION IN 200 CYCLES ***
        # Note that some lines have '(*******%)' in them, so we add spaces around
        # The asterisks to avoid triggering early.
        while " *** " not in line and "GEOMETRY OPTIMIZATION" not in line:
            if "* total energy  :" in line:
                line_split = line.split()
                energy_scf = float(line_split[4])
                self.parsed_info["outputs"]["energy_scf"].append(energy_scf)
            line = next(f)

    def md_avg_props(self, f, line):
        r"""Average MD properties.

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

            average properties
            Epot               :  -991.673089209627
            Epot (accurate SCC):  -991.771829504958
            Ekin               :   1.79432794141482
            Etot               :  -989.878761268212
            T                  :   522.456878068899
        """
        while "T                  :" not in line:
            if "Etot               :" in line:
                self.parsed_info["outputs"]["avg_energy_tot"] = float(line.split()[-1])
            line = next(f)

        self.parsed_info["outputs"]["avg_temp"] = float(line.split()[-1])

    def thermostat_prob(self, f, line):
        r"""If there was an issue with the thermostat.

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

            thermostating problem
        """
        self.parsed_info["runtime_info"]["success"] = False

    def success(self, f, line):
        r"""If the calculation is successful.

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

            normal termination of xtb
        """
        if "success" not in self.parsed_info["runtime_info"]:
            self.parsed_info["runtime_info"]["success"] = True
