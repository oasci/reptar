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


def prep_xtb_constrain_block(constraints):
    """Prepare xTB constrain block lines for input file.

    Parameters
    ----------
    constraints : :obj:`tuple`
        Internal constraints in a nested tuple to add to the xtb input file. The first
        element is the label (e.g., ``"distance"`` or ``"dihedral"``) and the second
        element is a list containing the values after  ``:`` used to control xtb.
        Atom indices should start from ``0``. For more information, see
        https://xtb-docs.readthedocs.io/en/latest/xcontrol.html#constraining-potentials

    Returns
    -------
    :obj:`list`
        Lines for the input file. Does not contain newline characters.
    """
    lines = []
    lines.append("$constrain")
    for constraint in constraints:
        constraint_name = constraint[0]
        if constraint_name == "force constant":
            constraint_line = f"    force constant={constraint[1]}"
        else:
            constraint_line = f"    {constraint[0]}: "
            # xtb indices start at 1
            constraint_line += ", ".join([str(i + 1) for i in constraint[1][:-1]])
            constraint_line += f", {constraint[1][-1]}"
        lines.append(constraint_line)
    lines.append("$end")
    return lines


def prep_xtb_input_lines(charge, multiplicity, constraints=None, save_traj=False):
    """Prepare lines for xtb input file.

    Parameters
    ----------
    charge : :obj:`int`
        Total charge of the system.
    multiplicity : :obj:`int`
        Spin state multiplicity of the system.
    constraints : :obj:`tuple`
        Internal constraints in a nested tuple to add to the xtb input file. The first
        element is the label (e.g., ``"distance"`` or ``"dihedral"``) and the second
        element is a list containing the values after  ``:`` used to control xtb.
        Atom indices should start from ``0``. For more information, see
        https://xtb-docs.readthedocs.io/en/latest/xcontrol.html#constraining-potentials

    Returns
    -------
    :obj:`list`
        Lines of xtb input file.
    """
    spin = int(multiplicity - 1)
    xtb_input_lines = [f"$chrg {charge} $end", f"$spin {spin} $end"]
    if constraints is not None:
        constraint_lines = prep_xtb_constrain_block(constraints)
        xtb_input_lines.extend(constraint_lines)
    if save_traj:
        xtb_input_lines.extend(["$opt", "    logfile=xtbopt.trj", "$end"])
    xtb_input_lines = [i + "\n" for i in xtb_input_lines]
    return xtb_input_lines
