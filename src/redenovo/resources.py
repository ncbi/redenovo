import os
import sys
import pandas as pd
import numpy as np

#get logger instance
import logging
_logger = logging.getLogger(__name__)


class Resources(object):

    def __init__(self, args, COSMIC):
        """ Initializes all tensors by either reading content from file, or
        creating random initial data for it.

        """
        _logger.debug("Creating Resources instance.")

        self.args = args
        self.COSMIC = COSMIC
        # We Store the tensors as named pandas dataframes internally so we can export them later
        # Numpy arrays do not natively support row names.

        # format for M
        # TODO add format description
        self.M = None

        # format for P
        # since P can be a mixture of fixed and unknown signatures,
        # we need to account for it by splitting it into these two
        # components. Hence {'fixed':[None or df], 'inferred':[None or df]}
        self.P = None
        self.A = None
        
        _logger.debug("Initializing tensors.")
        self.initialize_m()
        self.initialize_p()
        self.initialize_a()
        _logger.debug("Tensor initialization complete.")
        _logger.debug("Final tensor shapes are [dimensions with a \"+\" refer to the fixed (left) and the to be inferred component (right)]:")
        _logger.info(f"M: [{self.M.shape[0]}, {self.M.shape[1]}]")
        _logger.info(f"P: [{0 if self.P['fixed'] is None else self.P['fixed'].shape[0]} + {0 if self.P['inferred'] is None else self.P['inferred'].shape[0]},  {self.P['fixed'].shape[1] if self.P['fixed'] is not None else self.P['inferred'].shape[1]}]")
        _logger.info(f"A: [{self.A.shape[0]}, {self.A.shape[1]}]")

        # store performance values of interest here.
        # currently losses, best_losses, epochs, and steps
        self.logs = {}

    def initialize_m(self):
        """Depending on the parameters specified by the user, load M from file(s)
        """

        _logger.debug("Starting reading M")

        # First, make sure the file specified for M exist on non-volatile memory
        file = os.path.abspath(self.args.matrix)

        if not os.path.exists(file):
            _logger.critical(f"File '{os.path.basename(file)}' for input matrix M could not be found. Exiting.")
            sys.exit(1)

        _logger.debug(f"Parsing file {os.path.basename(file)}.")

        # Next, read in the file
        self.M = {}

        if self.args.has_header_and_index:
            m = pd.read_csv(
                file,
                sep=self.args.delimiter,
                header=0,      # First row as column headers
                index_col=0    # First column as row names
            )
            
            if set(list(self.COSMIC.columns)) != set(list(m.columns)):
                raise ValueError(f"Column (SBS) mismatch between COSMIC and input matrix M")

            m = m[self.COSMIC.columns]
            
        else:
            m = pd.read_csv(file,
                            sep=self.args.delimiter,
                            header=None,
                            index_col=False)
            m.columns = self.COSMIC.columns
            m.index = [f"genome_{x+1}" for x in range(m.shape[0])]

        # Cast to float
        m = m.astype(np.float64, copy=False, errors='raise')

        if self.args.bootstrap:
            m = m.sample(frac=0.8)#, random_state=1)

        self.M = m
        _logger.debug(f"Successfully read M with shape [{self.M.shape[0]},{self.M.shape[1]}]")

    def initialize_p(self):
        """Depending on the parameters specified by the user, load P from file(s),
        and/or initialize random matrix according to NUMPRI

        At the end of this fuction, either both or at least one of the options for
        Q will be present.

        Args:
          args: parser arguments
        """

        _logger.debug("Starting initizalize_p")

        # we need to enforce that M is read from file first, as the dimension
        # compatibility depends on it
        if self.M is None:
            _logger.debug("initialize_m has not been called yet. Make sure you call this function first.")
            _logger.critical("Internal error (no M). Exiting.")
            sys.exit(1)

        self.P = {}  # initialize

        if self.args.primary is not None:
            p = self.COSMIC.loc[self.args.primary]
            # cast to float
            p = p.astype(np.float64, copy=False, errors='raise')
            # Sanity check: the number of mutational categories must be identical to the one in M
            if self.M.shape[1] != p.shape[1]:
                _logger.critical(f"The number of mutational signatures of M and P must be identical. Currently, M has {self.M[0][1].shape[0]}, and P has {p.shape[1]} signatures. Exiting.")
                sys.exit(1)

            # add to dictionary
            self.P['fixed'] = p
            _logger.debug(f"Successfully read P with shape [{p.shape[0]},{p.shape[1]}]")
            
        else:  # no primary specified
            self.P['fixed'] = None

        # next, check if we have primary signatures to infer
        if self.args.numpri is not None:
            _logger.debug("Case N is not None")

            # initialize with random data
            p = pd.DataFrame(np.random.rand(self.args.numpri, self.M.shape[1]), # since the dimension of p in this case depends on knowing the number of categories, we source it from M
                                columns=self.M.columns,
                                index=[f"ReDeNovo_p{x+1}" for x in range(self.args.numpri)]
                            )

            # cast to float
            p = p.astype(np.float64, copy=False, errors='raise')

            # add to dictionary
            self.P['inferred'] = p

            _logger.debug(f"Successfully created {self.args.numpri} primary signatures to be infered with shape [{p.shape[0]},{p.shape[1]}]")

        else:  # no N
            self.P['inferred'] = None

        # Final report
        p_dim1 = (0 if self.P['fixed'] is None else self.P['fixed'].shape[0]) + (0 if self.P['inferred'] is None else self.P['inferred'].shape[0])
        p_dim2 = self.P['fixed'].shape[1] if self.P['fixed'] is not None else self.P['inferred'].shape[1]
        _logger.debug(f"Final shape for P is [{p_dim1},{p_dim2}]")

    def initialize_a(self):
        """Initializes A
        """

        _logger.debug("Starting initialize_a")

        # Sanity check. We need M and P to be present in order to determine the dimension of W
        if self.M is None or self.P is None:
            _logger.debug("initialize_m or initialize_p has not been called yet. Make sure you call this function first.")
            _logger.critical("Internal error (no M or P).  Exiting.")
            sys.exit(1)

        # determine the number of primary signatures and their row names from P
        num_pri_sigs = (0 if self.P['fixed'] is None else self.P['fixed'].shape[0]) + (0 if self.P['inferred'] is None else self.P['inferred'].shape[0])
        pri_sigs_index = ([] if self.P['fixed'] is None else list(self.P['fixed'].index)) + ([] if self.P['inferred'] is None else list(self.P['inferred'].index))

        # initialize with random data
        a = pd.DataFrame(   np.random.rand(len(self.M), num_pri_sigs), # since the dimension of p in this case depends on knowing the number of categories, we source it from M
                            columns=pri_sigs_index,
                            index=self.M.index)

        # cast to float
        a = a.astype(np.float64, copy=False, errors='raise')

        # add to dictionary
        self.A = a

        _logger.debug(f"Successfully created A with shape [{a.shape[0]},{a.shape[1]}]")

    def _get_dynamic_tensor(self, t, t_name="tensor"):
        """ Generic getter function for P
        Returns an array of numpy objects ['fixed', 'inferred'] with the content of t

            Args:
                t: The 2-dimensional tensor to return
                t_name: string representing the name of t for output purposes

            Returns:
                array with [ 'fixed' or None, 'inferred' or None]
        """
        _logger.debug(f"_get_dynamic_tensor() for {t_name}")

        # Sanity check. Has P been initialized?
        if t is None:
            _logger.critical(f"Tensor {t_name} has not yet been initialized. Call the corresponding initialize function first. Exiting")
            sys.exit(1)

        return [ None if t['fixed'] is None else t['fixed'].to_numpy(copy=True), None if t['inferred'] is None else t['inferred'].to_numpy(copy=True)]

    def _set_dynamic_tensor(self, T, t, t_name="tensor"):
        """ Generic setter for P.
        Sets the content of m to the corresponding pandas tables in M

        Args:
            T: Pandas dataframe to be set as t
            t: The 2-dimensional tensor to set
            t_name: string representing the name of t for output purposes

        """
        _logger.debug(f"_set_dynamic_tensor() for {t_name}")

        # Sanity check. Has M been initialized?
        if T is None:
            _logger.critical(f"Tensor {t_name} has not yet been initialized. Call initialize first. Exiting")
            sys.exit(1)

        # Sanity check. Are the dimensions of p the same as P?
        if T['fixed'] is None and t[0] is not None:
            _logger.critical(f"Trying to set {t_name} with a shape different from the one present.")
            _logger.debug(f"Passed shape for fixed component is {list(t[0].shape)} whereas stored shape is None. Exiting.")
            sys.exit(1)

        if T['fixed'] is not None and T['fixed'].shape != t[0].shape:
            _logger.critical(f"Trying to set {t_name} with a shape different from the one present.")
            _logger.debug(f"Passed shape for fixed component is {list(t[0].shape)} whereas stored shape is {T['fixed'].shape}. Exiting.")
            sys.exit(1)

        if T['inferred'] is None and t[1] is not None:
            _logger.critical(f"Trying to set {t_name} with a shape different from the one present.")
            _logger.debug(f"Passed shape for inferred component is {list(t[1].shape)} whereas stored shape is None. Exiting.")
            sys.exit(1)

        if T['inferred'] is not None and T['inferred'].shape != t[1].shape:
            _logger.critical(f"Trying to set {t_name} with a shape different from the one present.")
            _logger.debug(f"Passed shape for fixed component is {list(t[1].shape)} whereas stored shape is {T['inferred'].shape}. Exiting.")
            sys.exit(1)

        # Set values
        if t[0] is not None:
            T['fixed'][:] = t[0]

        if t[1] is not None:
            T['inferred'][:] = t[1]

    def _get_tensor(self, t, t_name="tensor"):
        """ Generic getter for A

            Args:
                t: The 2-dimensional tensor to return
                t_name: string representing the name of t for output purposes

            Returns:
                Copy of t as a numpy object
        """

        _logger.debug(f"getTensor() for {t_name}")

        # Sanity check. Has t been initialized?
        if t is None:
            _logger.critical(f"Tensor {t_name} has not yet been initialized. Call the appropriate initialize function first. Exiting")
            sys.exit(1)

        return t.to_numpy(copy=True)

    def _set_tensor(self, T, t, t_name="tensor"):
        """ Generic setter for A

            Args:
                T: Pandas dataframe to be set as t
                t: The 2-dimensional tensor to set
                t_name: string representing the name of t for output purposes

        """

        _logger.debug(f"setTensor() for {t_name}")

        # Sanity check. Has t been initialized?
        if t is None:
            _logger.critical(f"Tensor {t_name} has not yet been initialized. Call the appropriate initialize function first. Exiting")
            sys.exit(1)

        # Do shapes match?
        if T.shape != t.shape:
            _logger.critical(f"Trying to set {t_name} with a shape different from the one present.")
            _logger.debug(f"Passed shape is {list(t.shape)} whereas stored shape is {T.shape}. Exiting.")
            sys.exit(1)

        # Then set
        T[:] = t

    def get_m(self):
        """ Returns a numpy object with the content of M

            Returns:
                numpy with shape (GxK)
        """
        _logger.debug("get_m()")

        # Sanity check. Has M been initialized?
        if self.M is None:
            _logger.critical("M has not yet been initialized. Call initialize_m() first. Exiting")
            sys.exit(1)

        return self.M.to_numpy(copy=True)

    def set_m(self, m):
        """ Sets the content of m to the corresponding pandas tables in M
        """
        _logger.debug("set_m()")

        # Sanity check. Has M been initialized?
        if self.M is None:
            _logger.critical("M has not yet been initialized. Call initialize_m() first. Exiting")
            sys.exit(1)

        # Sanity check. Are the dimensions of m the same as M?
        if m.shape[0] != self.M.shape[0] or m.shape[1] != self.M.shape[1] or len(m.shape) != 2:
            _logger.critical("Trying to set M with a shape different from the one present.")
            _logger.debug(f"Passed shape is {list(m.shape)} whereas stored shape is [list({self.M.shape})]. Exiting.")
            sys.exit(1)

        # Set values
        for x in range(len(self.M)):
            self.M = m.copy()

    # define convecience functions for P, W, A, R, Q, and D.
    def get_p(self):
        return self._get_dynamic_tensor(self.P, 'P')

    def set_p(self, p):
        return self._set_dynamic_tensor(self.P, p, 'P')

    def set_p_fixed(self, p_fixed):
        p_cur = self.get_p()
        p_cur['fixed'] = p_fixed
        self.set_p(p_cur)
        return self._set_dynamic_tensor(self.P, p_cur, 'P with new fixed')

    def get_a(self):
        return self._get_tensor(self.A, 'A')

    def set_a(self, a):
        self._set_tensor(self.A, a, "A")

    def get_logs(self):
        return self.logs

    def set_logs(self, logs):
        for log in logs:
            self.logs[log] = logs[log]
