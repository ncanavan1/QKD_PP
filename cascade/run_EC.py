from cascade.algorithm import ALGORITHMS
from cascade.key import Key
from cascade.mock_classical_channel import MockClassicalChannel
from cascade.reconciliation import Reconciliation



def run_reconciliation(algorithm, key_size, error_method, error_rate,power_params):

    correct_key = Key.create_random_key(key_size)
    noisy_key = correct_key.copy(error_rate, error_method)
    actual_bit_errors = correct_key.difference(noisy_key)
    actual_bit_error_rate = actual_bit_errors / key_size
    mock_classical_channel = MockClassicalChannel(correct_key,power_params)
    reconciliation = Reconciliation(algorithm, mock_classical_channel, noisy_key, error_rate,power_params)
    reconciliated_key, Eves_traces = reconciliation.reconcile()
    remaining_bit_errors = correct_key.difference(reconciliated_key)

    return correct_key, noisy_key, reconciliated_key, Eves_traces, remaining_bit_errors