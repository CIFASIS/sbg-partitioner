import argparse
import json

def generate_reinits(var_id, conns, offset):
    records = []
    for neuron_conn in range(conns):
        record = {
            "id": var_id,
            "exp": [[1, offset+neuron_conn]],
            "defs": []
        }
        records.append(record)
    return records

def generate_neuron_connection_records(conns, size):
    offset = 100
    records = []
    # Generate list for partial total node definition
    for ev in range(2):
        record = {
            "id": int(7+ev),
            "interval": [[1, size-conns-offset]],
            "lhs": [
                {
                    "id": "ev_" + str(7+ev),
                    "exp": [[1, 0]],
                    "defs": []
                },
                {
                    "id": "v",
                    "exp": [[1, 0]],
                    "defs": []
                },
                {
                    "id": "ginh",
                    "exp": [[1, 0]],
                    "defs": []
                },
                {
                    "id": "gex",
                    "exp": [[1, 0]],
                    "defs": []
                },
                {
                    "id": "active",
                    "exp": [[1, 0]],
                    "defs": []
                },
                {
                    "id": "tfire",
                    "exp": [[1, 0]],
                    "defs": []
                }
            ],
            "rhs": [
                {
                    "id": "v",
                    "exp": [[1, 0]],
                    "defs": [1, 8, 9]
                }
            ]
        }
        if ev == 1:
            record["lhs"] += generate_reinits("exReinit", conns, offset)
        else:
            record["lhs"] += generate_reinits("inhReinit", conns, offset)
        records.append(record)

    return json.dumps(records)

def generate_json(size, conns, inputs):
    """
    Generates a JSON file input for the spiking neurons Modelica model.

    Parameters:
    size (int): The number of neurons.
    conns (int): The number of connections between neurons.
    inputs (int): The number of noise inputs in the network.

    Returns:
    str: A JSON-formatted string for spiking neurons Modelica model used as SBG Partitioner input.
    """
    data = {
        "nodes": [
            {
                "id": 1,
                "interval": [[1, size]],
                "lhs": [
                    {
                        "id": "v",
                        "exp": [[1, 0]],
                        "defs": []
                    }
                ],
                "rhs": [
                    {
                        "id": "active",
                        "exp": [[1, 0]],
                        "defs": [6, 8, 9]
                    },
                    {
                        "id": "v",
                        "exp": [[1, 0]],
                        "defs": [1, 8, 9]
                    },
                    {
                        "id": "gex",
                        "exp": [[1, 0]],
                        "defs": [2, 4, 8, 9]
                    },
                    {
                        "id": "ginh",
                        "exp": [[1, 0]],
                        "defs": [3, 5, 8, 9]
                    }
                ]
            },
            {
                "id": 2,
                "interval": [[1, size]],
                "lhs": [
                    {
                        "id": "gex",
                        "exp": [[1, 0]],
                        "defs": []
                    }
                ],
                "rhs": [
                    {
                        "id": "gex",
                        "exp": [[1, 0]],
                        "defs": [2, 4, 8, 9]
                    }
                ]
            },
            {
                "id": 3,
                "interval": [[1, size]],
                "lhs": [
                    {
                        "id": "ginh",
                        "exp": [[1, 0]],
                        "defs": []
                    }
                ],
                "rhs": [
                    {
                        "id": "ginh",
                        "exp": [[1, 0]],
                        "defs": [3, 5, 8, 9]
                    }
                ]
            },
            {
                "id": 4,
                "interval": [[1, size]],
                "lhs": [
                    {
                        "id": "ev_1",
                        "exp": [[1, 0]],
                        "defs": []
                    },
                    {
                        "id": "exReinit",
                        "exp": [[1, 0]],
                        "defs": []
                    },
                    {
                        "id": "gex",
                        "exp": [[1, 0]],
                        "defs": []
                    }
                ],
                "rhs": [
                    {
                        "id": "exReinit",
                        "exp": [[1, 0]],
                        "defs": [4, 7, 8]
                    },
                    {
                        "id": "gex",
                        "exp": [[1, 0]],
                        "defs": [2, 4, 8, 9]
                    }
                ]
            },
            {
                "id": 5,
                "interval": [[1, size]],
                "lhs": [
                    {
                        "id": "ev_2",
                        "exp": [[1, 0]],
                        "defs": []
                    },
                    {
                        "id": "inhReinit",
                        "exp": [[1, 0]],
                        "defs": []
                    },
                    {
                        "id": "ginh",
                        "exp": [[1, 0]],
                        "defs": []
                    }
                ],
                "rhs": [
                    {
                        "id": "inhReinit",
                        "exp": [[1, 0]],
                        "defs": [5, 8]
                    },
                    {
                        "id": "ginh",
                        "exp": [[1, 0]],
                        "defs": [3, 5, 8, 9]
                    }
                ]
            },
            {
                "id": 6,
                "interval": [[1, size]],
                "lhs": [
                    {
                        "id": "ev_3",
                        "exp": [[1, 0]],
                        "defs": []
                    },
                    {
                        "id": "active",
                        "exp": [[1, 0]],
                        "defs": []
                    }
                ],
                "rhs": [
                    {
                        "id": "tfire",
                        "exp": [[1, 0]],
                        "defs": [8, 9]
                    }
                ]
            },
            {
                "id": 7,
                "interval": [[1, inputs]],
                "lhs": [
                    {
                        "id": "ev_4",
                        "exp": [[1, 0]],
                        "defs": []
                    },
                    {
                        "id": "exReinit",
                        "exp": [[1, 100]],
                        "defs": []
                    },
                    {
                        "id": "inputs",
                        "exp": [[1, 0]],
                        "defs": []
                    },
                    {
                        "id": "nextev",
                        "exp": [[1, 0]],
                        "defs": []
                    }
                ],
                "rhs": [
                    {
                        "id": "inputs",
                        "exp": [[1, 0]],
                        "defs": [7]
                    },
                    {
                        "id": "nextev",
                        "exp": [[1, 0]],
                        "defs": [7]
                    }
                ]
            }
        ]
    }

    # Concatenate section controller records
    connections_records = json.loads(generate_neuron_connection_records(conns, size))
    data["nodes"].extend(connections_records)

    return json.dumps(data)

def main():
    parser = argparse.ArgumentParser(description='Generate SBG JSON input from spiking neurons Modelica model.')
    parser.add_argument('size', type=int, help='The number of neurons population.')
    parser.add_argument('conns', type=int, help='The number of connections between neurons.')
    parser.add_argument('inputs', type=int, help='The number of noise inputs in the network.')

    args = parser.parse_args()

    result = generate_json(args.size, args.conns, args.inputs)
    print(result)

if __name__ == '__main__':
    main()