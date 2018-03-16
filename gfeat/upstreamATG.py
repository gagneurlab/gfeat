class UpstreamATG:
    def __init__(self, allow_ORF = True, verbose_output = False):
        self.allow_ORF = allow_ORF
        self.verbose_output = verbose_output
        pass

    def predict_on_sample(self, seq):
        """
        :param seq: string utr's sequence
        :return: if verbose_output: dictionary:
                     first entry – 1 or 0 depending whether the uAUG is in-frame or not
                     second – 1 or 0 depending whether it corresponds to a uORF or not
                 else: NumPy array of 1 and 0 depending whether the uAUG is in-frame or not
        :example: if the input 5'UTR has 5 AUG, then
                    {
                    "frame": [1, 1, 0, 0, 1],
                    "uORF": [1, 1, 1, 0, 0]
                    }
        """
        if self.verbose_output:
            pass
        else:
            pass

    def predict_on_batch(self, seq_list):
        """
        :param seq_list: list of string utr's sequences
        :return: if verbose_output: NumPy array of dictionaries:
                    first entry – 1 or 0 depending whether the uAUG is in-frame or not
                    second – 1 or 0 depending whether it corresponds to a uORF or not
                 else: NumPy array of 1 and 0 whether the uAUG is in-frame or not
        """
        if self.verbose_output:
            pass
        else:
            pass
