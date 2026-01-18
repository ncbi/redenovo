import configargparse
import logging
import os
import pkg_resources
import json

from redenovo import __version__

__author__ = "ReDeNovo"
__copyright__ = "ReDeNovo"
__license__ = "GPL-3.0-only"

_logger = logging.getLogger(__name__)


def check_genome(value):
    """Validate integer input in arguments.

    Args:
      value: input value to check
      minval: value must be equal to or larger than minval

    Returns:
      True is sanity check passes, raises error otherwise
    """
    try:
        ivalue = int(value)
    except:
        raise configargparse.ArgumentTypeError(f"{value} is not an integer")

    if ivalue < 37 or ivalue > 38:
        raise configargparse.ArgumentTypeError(f"{value} is an invalid int value. Must be 37 or 38")

    return ivalue


def check_positive_int_or_zero(value, minval=0):
    """Validate integer input in arguments.

    Args:
      value: input value to check
      minval: value must be equal to or larger than minval

    Returns:
      True is sanity check passes, raises error otherwise
    """
    try:
        ivalue = int(value)
    except:
        raise configargparse.ArgumentTypeError(f"{value} is not an integer")

    if ivalue < minval:
        raise configargparse.ArgumentTypeError(f"{value} is an invalid int value. Must be >= {minval}")

    return ivalue


def check_positive_int_or_zero_or_minusone(value, minval=-1):
    """Validate integer input in arguments.

    Args:
      value: input value to check
      minval: value must be equal to or larger than minval

    Returns:
      True is sanity check passes, raises error otherwise
    """
    try:
        ivalue = int(value)
    except:
        raise configargparse.ArgumentTypeError(f"{value} is not an integer")

    if ivalue < minval:
        raise configargparse.ArgumentTypeError(f"{value} is an invalid int value. Must be >= {minval}")

    return ivalue


def check_positive_int(value, minval=1):
    """Validate integer input in arguments.

    Args:
      value: input value to check
      minval: value must be equal to or larger than minval

    Returns:
      True is sanity check passes, raises error otherwise
    """
    try:
        ivalue = int(value)
    except:
        raise configargparse.ArgumentTypeError(f"{value} is not an integer")

    if ivalue < minval:
        raise configargparse.ArgumentTypeError(f"{value} is an invalid int value. Must be >= {minval}")

    return ivalue

    
def check_positive_float(value, minval=0.0):
    """Validate integer input in arguments.

    Args:
      value: input value to check
      minval: value must be equal to or larger than minval

    Returns:
      True is sanity check passes, raises error otherwise
    """
    try:
        fvalue = float(value)
    except:
        raise configargparse.ArgumentTypeError(f"{value} is not a float")

    if fvalue < minval:
        raise configargparse.ArgumentTypeError(f"{value} is an invalid float value. Must be >= {minval}")

    return fvalue


def check_lower(value):
    """Converts the string to lower case

    Args:
        value: input value to convert

    Returns:
        Lower case represntation of value
    """

    return value.lower()


def check_upper(value):
    """Converts the string to upper case

    Args:
        value: input value to convert

    Returns:
        Upper case represntation of value
    """

    return value.upper()


def get_default_config_path():
    """Return the path on the file systems
    where the default config file is located as
    this might differ depending on where the
    package was installed.

    See: https://setuptools.readthedocs.io/en/latest/pkg_resources.html#basic-resource-access
    """

    resource_package = __name__
    resource_path = '/'.join(('data', 'default.conf'))
    file = pkg_resources.resource_filename(resource_package, resource_path)

    return file


def parse_args(args):
    """Parse command line parameters

    Args:
      args (List[str]): command line parameters as list of strings
          (for example  ``["--help"]``).

    Returns:
      :obj:`argparse.Namespace`: command line parameters namespace
    """
    
    parser = configargparse.ArgumentParser(default_config_files=[get_default_config_path()])

    #mandatory mutation count matrix
    parser.add('-M', '--matrix', required=True,  help='Mutation count matrix, specified as mutational profiles.')

    #primary signatures. either P or N or both are required, but at least one
    p_group = parser.add_argument_group('Primary signature options (at least one option required)')
    #p_group.add_argument('-P', '--primary', help='File containing the primary signature matrix of size NxK. If omitted, primary signatures will be infered according to -N.')
    p_group.add_argument('-P', '--primary', type=str, nargs='+', default=[], help='List of SBS values, e.g. ["SBS1", "SBS5"]. Default: []')
    p_group.add_argument('-N', '--numpri', type=check_positive_int_or_zero_or_minusone, default=-1, help='Number of novel signatures to infer. If -P is specificed in addition to -N, ReDeNovo will infer N signatures while keeping the ones defined in P as fixed. -1 for ReDenovo to decide it.')
    
    #additional input options
    #parser.add('-l', '--labels', default=False, action='store_true', help=f"Whether the input contains labels. If set, ReDeNovo will treat the first row/column of each input matrix as labels and use this information when generating the output files.")
    parser.add_argument('-g', '--genome',  type=check_genome, default=38, help='Genome version, either 37 or 38. Default=38.')
    parser.add_argument('-w', '--whole',  type=str, default='WGS', help='Sequencing platform for the data and COSMIC catalogue. "WGS" for Whole Genome Sequencing, "WES" for Whole Exome Sequencing. Default= "WGS"')
    parser.add_argument('--cosmic-version',  type=str, default='3.4', help='COSMIC version ["3.4", "3.3", "3.2", "3.1", "3", "2", or "1"]. Default= "3.4"')
    parser.add_argument('--manual-cosmic', action='store_true', help="Whether use manual COSMIC [user-provided COSMIC.txt in the input folder] (True or False). Default: False")
    parser.add_argument('--manual-cosmic-file', type=str, default=None, help="Path to the file containing reference signatures. Used only if --manual-cosmic is set.")
    parser.add_argument('--exposure-thr1', default=0.05, type=check_positive_float, help='Minimum patient-wise normalized exposure required for a signature to be considered present (default: 0.05)') 
    parser.add_argument('--exposure-thr2', default=1, type=check_positive_int, help='Minimum raw exposure required for a signature to be considered present (default: 1)')
    parser.add_argument('--exposure-thr3', default=100, type=check_positive_int, help='Minimum raw exposure required for a signature to be considered present for binary exposure matrix (default: 100)')
    parser.add_argument('--thr1', default=0.1, type=check_positive_float, help='Minimum fraction of patients with exposure ≥ thr1 required for a signature to be considered present (default: 0.1)')
    parser.add_argument('--thr2', default=0.70, type=check_positive_float, help='Minimum cosine similarity to match a signature with a known COSMIC signature and include in the inferred set (default: 0.70)')
    parser.add_argument('--thr3', default=0.70, type=check_positive_float, help='Minimum exposure weight for a signature to contribute to the final exposure profile (default: 0.70)') 
    parser.add_argument('--thr4', default=0.80, type=check_positive_float, help='Minimum cosine similarity to consider a signature as known and exclude it from novel candidate detection (default: 0.80)') 
    parser.add_argument('--thr5', default=0.1, type=check_positive_float, help='Minimum fraction of the cohort with nonzero exposure required for a signature to be considered present (default: 0.1)')
    parser.add_argument('-n', '--numruns', default=10, type=check_positive_int, help='Number of runs to repeat. Default=10') 
    parser.add_argument('-i', '--numiters', default=10, type=check_positive_int, help='Maximum number of iterations allowed while attempting to add new fixed signature (default: 10)') 
    parser.add_argument('--consno', default=1, type=check_positive_int, help='Minimum number of times a signature must be selected to be included in the inferred signature set (default: 1)') 
    p_group.add_argument('-E', '--exclude', type=str, nargs='+', default=[], help='List of SBS values to exclude from COSMIC, e.g. ["SBS1", "SBS5"]. Default: []')
    parser.add('-d', '--delimiter', default='\t', help="The delimiter used to separate the column in the input matrices. Default is tabulator. This delimiter will also be used for the output files.")

    #output options
    parser.add_argument('-O', '--out', required=False, default=os.getcwd(), help='Path to output folder. Folder will be created if it does not not exist. If omitted, the current folder will be used as the output directory. Existing files will be overwritten.')
    parser.add_argument('-v', '--verbosity', type=check_upper, default='INFO', choices=["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"], help='Set the verbosity level defining the detail of ReDeNovo messages to stdout.')
    parser.add_argument('--has-header-and-index', type=bool, default=False, help="Whether the file has row names and column headers (True or False)")
    parser.add_argument('--add-novel-signatures', action='store_true', help="Whether evaluate with novel signatures in the file provided (True or False). Used only if -N is set as 0. Default: False")
    parser.add_argument('--check-novel', action='store_true', help="Whether check the novel signature or not (True or False). Default: False")
    parser.add_argument('--novel-signatures-file', type=str, default=None, help="Path to the file containing novel signatures. Used only if --add-novel-signatures is set.")

    # optimizer options
    q_group = parser.add_argument_group('Optimizer options (all or none are required)')
    q_group.add_argument('-e', '--epochs', nargs='+', type=check_positive_int, help='List of integers specifying the number of epochs in each optimization bracket. List must be of same size as -s and -o')
    q_group.add_argument('-s', '--stepsizes', nargs='+', type=check_positive_float, help=f"List of integers specifying the stepsize for the corresponding bracket defined with -e. List must be of same size as -e and -o")
    q_group.add_argument('-o', '--optimizers', nargs='+', type=check_lower, choices=["adadelta", "adagrad", "adam", "adamax", "nadam", "rmaprop", "sgd"], help=f"List of optimizers to use in each bracket defined with -e. List must be of same size as -e and -s")

    # config file
    parser.add('-c', '--config', required=False, is_config_file=True, help='Config file path')

    # hidden arguments, not visiable to the user but used internally by redenovo
    parser.add_argument('--optimizer_user_update_steps', type=check_positive_int, help=configargparse.SUPPRESS)
    parser.add_argument('--optimizer_log_update_steps',  type=check_positive_int, help=configargparse.SUPPRESS)
    parser.add_argument('-m', '--misc', type=json.loads, help=configargparse.SUPPRESS)  # dictionary string for misc configs, '{"key":"value"}'
    #parser.add_argument('-r', '--run', type=check_positive_int, help=configargparse.SUPPRESS) 
    parser.add_argument('-b', '--bootstrap', default=False, help=configargparse.SUPPRESS) 
    #choices=['red', 'green', 'blue']

    args = parser.parse_args(args)

    # Additional sanity checks that are out of scope for the configargparse package
    # at least one primary signature option -P or -N
    if args.primary is None and args.numpri is None:
       parser.error("At least one primary signature option is required. Either use -P or -N, or both.")

    # N must be greater or equal to 1
    if args.numpri is not None and args.numpri < -1:
        parser.error("The number of primary signatures (-N) must be >= -1.")

    intersect_SBS = [i for i in args.exclude if i in args.primary]
    if len(intersect_SBS) > 0:
        parser.error("Excluded SBS should not exist in given primary signature set.")

    # Either all or none of the optimizer parameters must be specified. If they are, the lists must be of equal size
    opt_group = [args.epochs, args.stepsizes, args.optimizers]

    if sum(x is None for x in opt_group) not in [0, 3]:
        parser.error("Either all or non of the optimizer parameters (-e, -s, and -o) must be specified.")

    if sum(x is not None for x in opt_group) == 3 and len(set(map(len, opt_group))) != 1:
        parser.error(f"The optimizer parameters must specify the same number of brackets. Currently the number of brackets for each option is as follows: epochs:{len(args.epochs)}, stepsizes:{len(args.stepsizes)}, optimizers:{len(args.optimizers)}")

    return args
