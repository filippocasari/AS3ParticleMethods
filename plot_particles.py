import zmq
import numpy as np
import matplotlib.pyplot as plt

ctx = zmq.Context()
socket = ctx.socket(zmq.REP)
socket.bind("tcp://*:5555")

while True:
    # Receive first message
    message1 = socket.recv()
    array1 = np.frombuffer(message1, dtype=np.float64)

    # Receive second message
    message2 = socket.recv()
    array2 = np.frombuffer(message2, dtype=np.float64)

    # Plot arrays
    plt.scatter(array1, array2)
    plt.legend()
    plt.draw()

    # Send reply
    socket.send(b"OK")
