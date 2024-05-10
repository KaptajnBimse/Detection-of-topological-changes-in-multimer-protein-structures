import numpy as np
def find_increasing_subarrays(arr):
        # Initialize the current length and the result list
        current_length = 1
        result = []
        result2 = []

        # Iterate over the array
        for i in range(1, len(arr)):
            # If the current number is one greater than the previous number, increase the current length
            if arr[i] == arr[i - 1] + 1:
                current_length += 1
            else:
                # Otherwise, add the current length to the result list current_length times, and reset it
                result.extend(np.linspace(1, current_length, current_length, dtype=int))
                result2.extend([current_length]*current_length)
                current_length = 1

        # Don't forget to add the last subarray length
        result.extend(np.linspace(1, current_length, current_length, dtype=int))
        result2.extend([current_length]*current_length)

        return result, result2