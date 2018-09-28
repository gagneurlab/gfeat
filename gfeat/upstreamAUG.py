import re
import math
import numpy as np


class UpstreamAUG:
    def __init__(self, allow_ORF=True, verbose_output=False):
        """
        Constructor

        :param allow_ORF: bool, True by default, whether to check uORFs
        :param verbose_output: bool, False by default, whether to return dictionaries in predict_on_sample() and predict_on_batch() methods or not
        """
        self.allow_ORF = allow_ORF
        self.verbose_output = verbose_output
        pass

    def predict_on_sample(self, seq):
        """
        Predict_on_sample

        :param seq: string, 5'UTR's sequence

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

        if self.allow_ORF:
            if self.verbose_output:

                ATG_frame = []
                ATG_ORF = []

                for ATG in re.finditer('ATG', seq.upper()):
                    seq_remainder = seq[ATG.start() + 3:]
                    TAA_frame = [(TAA.start() % 3) for TAA in re.finditer('TAA', seq_remainder)]
                    if 0 in TAA_frame:
                        ORF = True
                    else:
                        TAG_frame = [(TAG.start() % 3) for TAG in re.finditer('TAG', seq_remainder)]
                        if 0 in TAG_frame:
                            ORF = True
                        else:
                            TGA_frame = [(TGA.start() % 3) for TGA in re.finditer('TGA', seq_remainder)]
                            ORF = 0 in TGA_frame
                    if ORF:
                        ATG_ORF.append(1)
                    else:
                        ATG_ORF.append(0)

                    if (len(seq) - ATG.start()) % 3:
                        ATG_frame.append(0)
                    else:
                        ATG_frame.append(1)
                return {"frame": np.array(ATG_frame), "uORF": np.array(ATG_ORF)}

            else:

                ATG_pos = [ATG.start() for ATG in re.finditer('ATG', seq.upper())]
                ATG_frame = [((len(seq) - pos) % 3) for pos in ATG_pos]
                ATG_frame[:] = [(math.ceil(res / 2) ^ 1) for res in ATG_frame]
                return np.array(ATG_frame)

        else:
            pass

    def predict_on_sample_with_pos(self, seq):
        """
        In comparison to predict_on_sample(), additionally returns the positions of AUGs

        :param seq: string utr's sequence

        :return: if verbose_output: dictionary

                first entry – 1 or 0 depending whether the uAUG is in-frame or not

                second – 1 or 0 depending whether it corresponds to a uORF or not

                third - pos of the ATG

            else: NumPy array of 1 and 0 depending whether the uAUG is in-frame or not

        :example: if the input 5'UTR has 5 AUG, then

            {

                "frame": [1, 1, 0, 0, 1],

                "uORF": [1, 1, 1, 0, 0],

                "pos": [38, 190, 438, 769, 981]

            }
        """

        if self.allow_ORF:
            if self.verbose_output:

                ATG_frame = []
                ATG_ORF = []
                ATG_pos = []

                for ATG in re.finditer('ATG', seq.upper()):
                    seq_remainder = seq[ATG.start() + 3:]
                    TAA_frame = [(TAA.start() % 3) for TAA in re.finditer('TAA', seq_remainder)]
                    if 0 in TAA_frame:
                        ORF = True
                    else:
                        TAG_frame = [(TAG.start() % 3) for TAG in re.finditer('TAG', seq_remainder)]
                        if 0 in TAG_frame:
                            ORF = True
                        else:
                            TGA_frame = [(TGA.start() % 3) for TGA in re.finditer('TGA', seq_remainder)]
                            ORF = 0 in TGA_frame
                    if ORF:
                        ATG_ORF.append(1)
                    else:
                        ATG_ORF.append(0)

                    if (len(seq) - ATG.start()) % 3:
                        ATG_frame.append(0)
                    else:
                        ATG_frame.append(1)
                    ATG_pos.append(ATG.start())
                return {"frame": np.array(ATG_frame), "uORF": np.array(ATG_ORF), "pos": np.array(ATG_pos)}

            else:

                ATG_pos = [ATG.start() for ATG in re.finditer('ATG', seq.upper())]
                ATG_frame = [((len(seq) - pos) % 3) for pos in ATG_pos]
                ATG_frame[:] = [(math.ceil(res / 2) ^ 1) for res in ATG_frame]
                return np.array(ATG_frame)

        else:
            pass

    def predict_on_sample_with_pos_pandas(self, seq, result_dict, strand, start=None):
        """
        In comparison to predict_on_sample(), additionally returns as positions of AUGs and outputs everything to the \
        passed to it dictionary

        :param seq: string utr's sequence
        :param result_dict: dictionary with 4 mandatory keys "not_in-frame_no_uORF", "not_in-frame_uORF", "in-frame_no_uORF", "in-frame_uORF", where to append the found values
        :param start: integer, position relatively to the whole genome (in contrast to position relative to the exon)

        """
        if self.allow_ORF:
            if strand == '+':
                if self.verbose_output:

                    list_00 = []  # not_in-frame_no_uORF
                    list_01 = []  # not_in-frame_uORF
                    list_10 = []  # in-frame_no_uORF
                    list_11 = []  # in-frame_uORF

                    for ATG in re.finditer('ATG', seq):
                        seq_remainder = seq[ATG.start() + 3:]
                        TAA_frame = [(TAA.start() % 3) for TAA in re.finditer('TAA', seq_remainder)]
                        if 0 in TAA_frame:
                            ORF = True
                        else:
                            TAG_frame = [(TAG.start() % 3) for TAG in re.finditer('TAG', seq_remainder)]
                            if 0 in TAG_frame:
                                ORF = True
                            else:
                                TGA_frame = [(TGA.start() % 3) for TGA in re.finditer('TGA', seq_remainder)]
                                ORF = 0 in TGA_frame
                        if ORF:
                            if (len(seq) - ATG.start()) % 3:
                                list_01.append(ATG.start() + start)
                            else:
                                list_11.append(ATG.start() + start)
                        else:
                            if (len(seq) - ATG.start()) % 3:
                                list_00.append(ATG.start() + start)
                            else:
                                list_10.append(ATG.start() + start)

                    result_dict["not_in-frame_no_uORF"].append(np.array(list_00))
                    result_dict["not_in-frame_uORF"].append(np.array(list_01))
                    result_dict["in-frame_no_uORF"].append(np.array(list_10))
                    result_dict["in-frame_uORF"].append(np.array(list_11))

                    pass

                else:

                    ATG_pos = [ATG.start() for ATG in re.finditer('ATG', seq)]
                    ATG_frame = [((len(seq) - pos) % 3) for pos in ATG_pos]
                    ATG_frame[:] = [(math.ceil(res / 2) ^ 1) for res in ATG_frame]
                    pass
            else:
                if self.verbose_output:

                    list_00 = []  # not_in-frame_no_uORF
                    list_01 = []  # not_in-frame_uORF
                    list_10 = []  # in-frame_no_uORF
                    list_11 = []  # in-frame_uORF

                    for ATG in re.finditer('ATG', seq):
                        seq_remainder = seq[ATG.start() + 3:]
                        TAA_frame = [(TAA.start() % 3) for TAA in re.finditer('TAA', seq_remainder)]
                        if 0 in TAA_frame:
                            ORF = True
                        else:
                            TAG_frame = [(TAG.start() % 3) for TAG in re.finditer('TAG', seq_remainder)]
                            if 0 in TAG_frame:
                                ORF = True
                            else:
                                TGA_frame = [(TGA.start() % 3) for TGA in re.finditer('TGA', seq_remainder)]
                                ORF = 0 in TGA_frame
                        if ORF:
                            if (len(seq) - ATG.start()) % 3:
                                list_01.append(start + (len(seq) - ATG.start()) - 1)
                            else:
                                list_11.append(start + (len(seq) - ATG.start()) - 1)
                        else:
                            if (len(seq) - ATG.start()) % 3:
                                list_00.append(start + (len(seq) - ATG.start()) - 1)
                            else:
                                list_10.append(start + (len(seq) - ATG.start()) - 1)

                    result_dict["not_in-frame_no_uORF"].append(np.array(list_00))
                    result_dict["not_in-frame_uORF"].append(np.array(list_01))
                    result_dict["in-frame_no_uORF"].append(np.array(list_10))
                    result_dict["in-frame_uORF"].append(np.array(list_11))

                    pass

                else:

                    ATG_pos = [ATG.start() for ATG in re.finditer('ATG', seq)]
                    ATG_frame = [((len(seq) - pos) % 3) for pos in ATG_pos]
                    ATG_frame[:] = [(math.ceil(res / 2) ^ 1) for res in ATG_frame]
                    pass
        else:
            pass

    def predict_on_sample_with_stop_pandas(self, seq, result_dict, strand, start=None):
        """
        In comparison to predict_on_sample(), additionally returns as positions of AUGs and outputs everything to the \
        passed to it dictionary

        :param seq: string utr's sequence
        :param result_dict: dictionary with 4 mandatory keys "not_in-frame_no_uORF", "not_in-frame_uORF", \
                            "in-frame_no_uORF", "in-frame_uORF", where to append the found values
        :param start: integer, position relatively to the whole genome (in contrast to position relative to the exon)

        """

        if self.allow_ORF:
            if strand == '+':
                if self.verbose_output:

                    list_00 = []  # not_in-frame_no_uORF
                    list_01 = []  # not_in-frame_uORF
                    list_10 = []  # in-frame_no_uORF
                    list_11 = []  # in-frame_uORF

                    for ATG in re.finditer('ATG', seq):
                        ORF = 0
                        seq_remainder = seq[ATG.start() + 3:]

                        for TAA in re.finditer('TAA', seq_remainder):
                            if not (TAA.start() % 3):
                                ORF = TAA.start()
                                break
                        if not ORF:
                            for TAG in re.finditer('TAG', seq_remainder):
                                if not (TAG.start() % 3):
                                    ORF = TAG.start()
                                    break
                            if not ORF:
                                for TGA in re.finditer('TGA', seq_remainder):
                                    if not (TGA.start() % 3):
                                        ORF = TGA.start()
                                        break
                        if ORF:
                            if (len(seq) - ATG.start()) % 3:
                                list_01.append(ATG.start() + start)
                                list_01.append(ORF + start)
                            else:
                                list_11.append(ATG.start() + start)
                                list_11.append(ORF + start)
                        else:
                            if (len(seq) - ATG.start()) % 3:
                                list_00.append(ATG.start() + start)
                            else:
                                list_10.append(ATG.start() + start)

                    result_dict["not_in-frame_no_uORF"].append(np.array(list_00))
                    result_dict["not_in-frame_uORF"].append(np.array(list_01))
                    result_dict["in-frame_no_uORF"].append(np.array(list_10))
                    result_dict["in-frame_uORF"].append(np.array(list_11))

                    pass

                else:

                    ATG_pos = [ATG.start() for ATG in re.finditer('ATG', seq)]
                    ATG_frame = [((len(seq) - pos) % 3) for pos in ATG_pos]
                    ATG_frame[:] = [(math.ceil(res / 2) ^ 1) for res in ATG_frame]
                    pass
            else:
                if self.verbose_output:

                    list_00 = []  # not_in-frame_no_uORF
                    list_01 = []  # not_in-frame_uORF
                    list_10 = []  # in-frame_no_uORF
                    list_11 = []  # in-frame_uORF

                    for ATG in re.finditer('ATG', seq):
                        ORF = 0
                        seq_remainder = seq[ATG.start() + 3:]

                        for TAA in re.finditer('TAA', seq_remainder):
                            if not (TAA.start() % 3):
                                ORF = TAA.start()
                                break
                        if not ORF:
                            for TAG in re.finditer('TAG', seq_remainder):
                                if not (TAG.start() % 3):
                                    ORF = TAG.start()
                                    break
                            if not ORF:
                                for TGA in re.finditer('TGA', seq_remainder):
                                    if not (TGA.start() % 3):
                                        ORF = TGA.start()
                                        break
                        if ORF:
                            if (len(seq) - ATG.start()) % 3:
                                list_01.append(start + (len(seq) - ATG.start()) - 1)
                                list_01.append(start + (len(seq) - ORF) - 1)
                            else:
                                list_11.append(start + (len(seq) - ATG.start()) - 1)
                                list_11.append(start + (len(seq) - ORF) - 1)
                        else:
                            if (len(seq) - ATG.start()) % 3:
                                list_00.append(start + (len(seq) - ATG.start()) - 1)
                            else:
                                list_10.append(start + (len(seq) - ATG.start()) - 1)

                    result_dict["not_in-frame_no_uORF"].append(np.array(list_00))
                    result_dict["not_in-frame_uORF"].append(np.array(list_01))
                    result_dict["in-frame_no_uORF"].append(np.array(list_10))
                    result_dict["in-frame_uORF"].append(np.array(list_11))

                    pass

                else:

                    ATG_pos = [ATG.start() for ATG in re.finditer('ATG', seq)]
                    ATG_frame = [((len(seq) - pos) % 3) for pos in ATG_pos]
                    ATG_frame[:] = [(math.ceil(res / 2) ^ 1) for res in ATG_frame]
                    pass
        else:
            pass

    def predict_on_batch(self, seq_list):
        """
        Predict on batch
        :param seq_list: list of string utr's sequences
        :return: if verbose_output: NumPy array of dictionaries

            first entry – 1 or 0 depending whether the uAUG is in-frame or not

            second – 1 or 0 depending whether it corresponds to a uORF or not

        else: NumPy array of 1 and 0 whether the uAUG is in-frame or not
        """

        if self.allow_ORF:

            result_list = []

            for seq in seq_list:
                result_list.append(self.predict_on_sample(seq))
            return result_list

        else:
            pass
