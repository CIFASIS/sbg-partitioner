import argparse

def transform_input_to_dict(input_string):
    # Remove unnecessary characters and split the string into partitions
    partitions = input_string.split('\\n')
    result_dict = {}
    for i, part in enumerate(partitions):
        part = part.strip().replace('{', '').replace('}', '')
        intervals = part.split('},')
        
        for interval in intervals:
            ranges = interval.strip().split(',')
            inner_list = [
                [int(start), int(end)] for r in ranges if r.strip() 
                for start, end in [map(int, r.strip('[]').split(':'))]
            ]
            result_dict[str(i)] = inner_list
    
    return result_dict


def get_next_index(input_list, current_index):
    # Return the next index or the last index if current is out of bounds
    return min(current_index + 1, len(input_list) - 1)


def has_valid_index(intervals, current_index):
    return current_index < len(intervals) - 1


def get_next_interval(result_dict, indexes, last_indexes):
    keys = list(result_dict.keys())
    
    current_min = float('inf')
    current_part = ""
    min_index = -1
    
    for i in range(len(keys)):
        part = keys[i]
        if not last_indexes[int(part)]:
            part_min = result_dict[part][indexes[i]][0]  # Get the start of the current range
            if part_min < current_min:
                min_index = indexes[i]
                current_min = part_min
                current_part = part
                
    return current_part, min_index


def is_valid_input(last_indexes):
    return not all(last_indexes)


def main():
    # Set up command line argument parsing
    parser = argparse.ArgumentParser(description='Process input string of intervals.')
    parser.add_argument('input_string', type=str, help='Input string of intervals in the format "{[start:end],[start:end]}"')
    
    args = parser.parse_args()
    input_data = args.input_string

    result = transform_input_to_dict(input_data)

    # Initialize indexes and last_indexes
    indexes = [0] * len(result)
    last_indexes = [False] * len(result)

    while is_valid_input(last_indexes):
        part, min_index = get_next_interval(result, indexes, last_indexes)
        
        begin, end = result[part][min_index]
        for _ in range(end - begin + 1):
            print(part)
        
        last_indexes[int(part)] = (len(result[part]) - 1 == min_index)
        indexes[int(part)] = get_next_index(result[part], indexes[int(part)])


if __name__ == "__main__":
    main()