from cascade.classical_channel import ClassicalChannel
import numpy as np

class MockClassicalChannel(ClassicalChannel):
    """
    A mock concrete implementation of the ClassicalChannel base class, which is used for the
    experiments.
    """

    def __init__(self, correct_key, power_params):
        self._correct_key = correct_key
        self._id_to_shuffle = {}
        self._reconciliation_started = False
        self._power_params = power_params

    def start_reconciliation(self):
        self._reconciliation_started = True

    def end_reconciliation(self):
        self._reconciliation_started = False
        self._id_to_shuffle = {}

    def ask_parities(self, blocks, Eves_traces):

        ##Eves information structured like 
        # Power, Bits, Outparity
        parities = []
        for block in blocks:
            
            shuffle = block.get_shuffle()
            start_index = block.get_start_index()
            end_index = block.get_end_index()
            parity, pwr_trace = shuffle.calculate_parity(self._correct_key, start_index, end_index,self._power_params)
            parities.append(parity)
            bit_pos_eve = [*shuffle._shuffle_index_to_key_index.values()][start_index:end_index]

           # if end_index - start_index < 2:
            #    k=7
            Eves_traces.append([pwr_trace,np.asarray(bit_pos_eve),parity])

        



        return parities, Eves_traces
