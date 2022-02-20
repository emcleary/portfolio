from typing import Callable, TypeVar, Union
import numpy.typing as npt
import numpy as np

NDArrayFloat32 = npt.NDArray[np.float32]
NDArrayFloat64 = npt.NDArray[np.float64]
NDArrayFloat = TypeVar('NDArrayFloat', NDArrayFloat32, NDArrayFloat64)

CallModel = Callable[[Union[NDArrayFloat32, NDArrayFloat64]], NDArrayFloat]

