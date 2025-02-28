from wisp.wisp_light.training.utils import softmax, softmax2
import numpy as np

def test_softmax():
    normalisation_func= "delta_mean"
    read_identity_threshold = 0.8
    predictions =np.ones(100, dtype=np.float32)
    result = softmax(predictions, normalisation_func, read_identity_threshold)
    print(result)
    result2 = softmax(predictions, normalisation_func, read_identity_threshold)

    print(result2)


test_softmax()