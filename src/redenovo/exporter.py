import os
# get logger instance
import logging
_logger = logging.getLogger(__name__)


class Exporter(object):

    def __init__(self, data):

        _logger.debug("Creating instance of Optimizer")

        self.data = data
        self.args = data.args
        self.outfolder = os.path.abspath(self.args.out)

        _logger.debug(f"Output folder is {self.outfolder}")

        os.makedirs(self.outfolder, exist_ok=True)
        _logger.debug(f"Output folder {self.outfolder} is ready.")
        # take care of the output folder
        #try:
        #    os.makedirs(self.outfolder)
        #except OSError:
        #    _logger.debug(f"Output folder {self.outfolder} exists.")

    def write_tables(self):
        """
        Writes the optimized tensors to file
        """
        _logger.info("Writing matrices to file")

        if self.data.P['fixed'] is not None:
            self.data.P['fixed'].to_csv(os.path.join(self.outfolder, 'P_fixed.txt'),
                                        sep=self.args.delimiter,
                                        index=True,
                                        header=True)

        if self.data.P['inferred'] is not None:
            self.data.P['inferred'].to_csv(os.path.join(self.outfolder, 'P_inferred.txt'),
                                            sep=self.args.delimiter,
                                            index=True,
                                            header=True)

        self.data.A.to_csv(os.path.join(self.outfolder, 'A.txt'),
                                        sep=self.args.delimiter,
                                        index=True,
                                        header=True)

    def write_tables_exposure(self):
        """
        Writes the exposure matrix to file
        """
        _logger.info("Writing matrices to file")

        A_binary = self.data.A[self.data.A < self.args.exposure_thr3] = 0
        A_binary.to_csv(os.path.join(self.outfolder, f'A_binary_with_thr_{self.args.exposure_thr3}.txt'),
                                        sep=self.args.delimiter,
                                        index=True,
                                        header=True)          
                                        
    # write optimization logs to file, usefull for plotting etc
    def write_logs(self):
        with open(os.path.join(self.outfolder, "optimization.log"), 'w') as flog:
            for log in self.data.logs:
                flog.write(f'{log}:{",".join(str(x) for x in self.data.logs[log])}\n')
