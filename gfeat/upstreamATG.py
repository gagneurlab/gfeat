# from gfeat.units import mutate_sequence
import re
import math
import numpy as np


class UpstreamATG:
    def __init__(self, allow_ORF=True, verbose_output=False):
        self.allow_ORF = allow_ORF
        self.verbose_output = verbose_output
        pass

    def predict_on_sample(self, seq):
        """
        NOTE: parameters have changed

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

        # mutated_seq = mutate_sequence(UTR_seq)
        #
        # mutated_exon_seq = []
        # ATG_pos = []
        # ATG_frame = []
        #
        # for seq in mutated_seq:
        #     temp_seq = ""
        #     temp_no = []
        #     for interval in exon_interval:
        #         temp_seq = temp_seq + seq[0][interval.start:interval.end]
        #         # Todo: decide what to do with numbers
        #     mutated_exon_seq.append(temp_seq)
        # for seq in mutated_exon_seq:
        #     ATG_pos.append([ATG.start() for ATG in re.finditer('ATG', seq.upper())])
        #     temp = [((len(seq) - pos) % 3) for pos in ATG_pos[-1]]
        #     temp[:] = [math.ceil(res / 2) for res in temp]
        #     ATG_frame.append(temp)

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
            pass  # Todo what should it do in this case?

    def predict_on_sample_with_pos(self, seq):
        """
        NOTE: parameters have changed

        :param seq: string utr's sequence
        :return: if verbose_output: dictionary:
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
            pass  # Todo what should it do in this case?

    def predict_on_sample_with_pos_pandas(self, seq, dict):
        """
        NOTE: parameters have changed

        :param seq: string utr's sequence
        :return: if verbose_output: dictionary:
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

                list_00 = []
                list_01 = []
                list_10 = []
                list_11 = []

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
                            list_01.append(ATG.start())
                        else:
                            list_11.append(ATG.start())
                    else:
                        if (len(seq) - ATG.start()) % 3:
                            list_00.append(ATG.start())
                        else:
                            list_10.append(ATG.start())

                dict["Transcript"].append(transcript_id)
                dict["0-0"].append(np.array(list_00))
                dict["0-1"].append(np.array(list_01))
                dict["1-0"].append(np.array(list_10))
                dict["1-1"].append(np.array(list_11))

                pass

            else:

                ATG_pos = [ATG.start() for ATG in re.finditer('ATG', seq)]
                ATG_frame = [((len(seq) - pos) % 3) for pos in ATG_pos]
                ATG_frame[:] = [(math.ceil(res / 2) ^ 1) for res in ATG_frame]
                return np.array(ATG_frame)

        else:
            pass  # Todo what should it do in this case?

    def predict_on_batch(self, seq_list):
        """
        :param seq_list: list of string utr's sequences
        :return: if verbose_output: NumPy array of dictionaries:
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
            pass  # Todo what should it do in this case?
