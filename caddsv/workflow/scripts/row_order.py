def inverse_permutation(order, expected_length):
    """Return source row positions in original input order."""
    if len(order) != expected_length:
        raise ValueError(
            f"Flank order has {len(order)} rows; expected {expected_length}."
        )
    if sorted(order) != list(range(expected_length)):
        raise ValueError("Flank order is not a permutation of the input rows.")

    inverse = [0] * expected_length
    for sorted_position, original_position in enumerate(order):
        inverse[original_position] = sorted_position
    return inverse


def read_inverse_permutation(path, expected_length):
    with open(path, encoding="utf-8") as order_file:
        order = [int(line.strip()) for line in order_file if line.strip()]
    return inverse_permutation(order, expected_length)
