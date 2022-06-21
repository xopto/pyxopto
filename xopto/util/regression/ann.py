# -*- coding: utf-8 -*-
################################ Begin license #################################
# Copyright (C) Laboratory of Imaging technologies,
#               Faculty of Electrical Engineering,
#               University of Ljubljana.
#
# This file is part of PyXOpto.
#
# PyXOpto is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PyXOpto is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PyXOpto. If not, see <https://www.gnu.org/licenses/>.
################################# End license ##################################

from typing import Tuple, List, Callable
import os.path
import json
import sys
import time

import numpy as np
import tensorflow.keras.models
import tensorflow.keras.layers
import tensorflow.keras.callbacks
import  tensorflow.keras.optimizers

from xopto.util.regression.pp import Preprocessor

TOPOLOGY_EXTENSION = '.json'
MODEL_EXTENSION = '.h5'


class Errors:
    @staticmethod
    def rmse(true: np.ndarray, estimated: np.ndarray) -> np.ndarray:
        '''
        Root mean square error estimator.

        Parameters
        ----------
        true: np.ndarray
            True values stored row-wise.
        estimated: np.ndarray
            Estimated values stored row-wise.

        Returns
        -------
        rmse: np.ndarray
            Root mean square error.
        '''
        true = np.asarray(true)
        if true.ndim == 1:
            true = np.expand_dims(true, -1)
        estimated = np.asarray(estimated)
        if estimated.ndim == 1:
            estimated = np.expand_dims(estimated, -1)

        return np.sqrt(((true - estimated)**2).mean(axis=0))

    @staticmethod
    def relative_rmse(true: np.ndarray, estimated: np.ndarray,
                      threshold: float = None) -> np.ndarray:
        '''
        Relative root mean square error estimator.

        Parameters
        ----------
        true: np.ndarray
            True values stored row-wise.
        estimated: np.ndarray
            Estimated values stored row-wise.
        threshold: float or None
            Only the true values that exceed this threshold are included
            in the evaluation.

        Returns
        -------
        rrmse: np.ndarray
            Estimated relative root mean square error (RMSE).
        '''
        true = np.asarray(true)
        if true.ndim == 1:
            true = np.expand_dims(true, -1)
        estimated = np.asarray(estimated)
        if estimated.ndim == 1:
            estimated = np.expand_dims(estimated, -1)

        if threshold is not None:
            ind = np.nonzero(true > threshold)
            true = true[ind]
            estimated = estimated[ind]

        return np.sqrt((((true - estimated)/true)**2).mean(axis=0))


class Topology:
    def __init__(self, inputs: int = 5,
                 node: int or Tuple[int] = None,
                 bias: bool or Tuple[bool] = True,
                 activation: str or Tuple[str] = 'relu',
                 initializer: str = 'normal',
                 pp_input: Preprocessor = None,
                 pp_output: Preprocessor = None):
        '''
        Sequential artificial neural network topology constructor.

        Parameters
        ----------
        inputs: int
            Number of inputs to the network (not the number of neurons in
            the first layer).
        node: int Tuple[int]
            Number of nodes/neurons in the layers. The number of elements in the
            list defines the number of layers. The last entry in the list
            represents the number of outputs, i.e. neurons in the output layer.
            The number of hidden layers is one less the number of elements in
            the neurons list/tuple.
        bias: bool or Tuple[bool]
            Enables additional bias for the layers. If a list or tuple of bools,
            the number of elements must match the number of layers. If a bool,
            the value applies to all the layers.
        activation: str or Tuple[str]
            Activation function of each layer. If a list or tuple of str,
            the number of elements must match the number of layers. If a str,
            the same activation function is used for all layers.
        initializer: str or Tuple[str]
            Generator for initialization of network parameters.
        pp_input: Preprocessor
            Preprocessor of the input variables.
        pp_output: Preprocessor
            Preprocessor of the output variables.
        '''
        self._inputs = int(inputs)

        self._model = None

        if isinstance(node, int):
            node = [node]
        if isinstance(node, (tuple, list)):
            node = list(node)
        for index, value in enumerate(node):
            node[index] = int(value)
        self._node = node

        n_layers = len(self._node)

        if not isinstance(activation, (list, tuple)):
            self._activation = [activation]*n_layers
        else:
            self._activation = activation

        if isinstance(bias, bool):
            self._bias = [bias]*n_layers
        else:
            self._bias = bias

        if not isinstance(initializer, (list, tuple)):
            self._initializer = [initializer]*n_layers
        else:
            self._initializer = initializer

        if pp_input is None:
            pp_input = Preprocessor()

        if pp_output is None:
            pp_output = Preprocessor()

        self._pp_input = pp_input
        self._pp_output = pp_output

    def _get_inputs(self) -> int:
        return self._inputs
    inputs = property(_get_inputs, None, None, 'Number of inputs.')

    def _get_node(self) -> Tuple[int]:
        return tuple(self._node)
    node = property(_get_node, None, None,
                    'The number of nodes/neurons used in the individual '
                    'layers.')

    def _get_bias(self) -> Tuple[str]:
        return tuple(self._bias)
    bias = property(_get_bias, None, None,
                    'The bias terms used in the individual layers.')

    def _get_activation(self) -> Tuple[str]:
        return self._activation
    activation = property(_get_activation, None, None,
                          'The activation functions used in the individual '
                          'layers.')

    def _get_initializer(self) -> Tuple[str]:
        return tuple(self._initializer)
    initializer = property(_get_initializer, None, None,
                           'The weight initialization method used in the '
                           'individual layers.')

    def _get_pp_input(self) -> Preprocessor:
        return self._pp_input
    pp_input = property(_get_pp_input, None, None,
                        'The input data preprocessor.')

    def _get_pp_output(self) -> Preprocessor:
        return self._pp_output
    pp_output = property(_get_pp_output, None, None,
                         'The output data preprocessor.')

    def model(self, verbose: bool = False) \
            -> tensorflow.keras.models.Sequential:
        '''
        Create an ANN model based on the specified topology.
        The last created model can be accessed through the public
        member .model.

        Parameters
        ----------
        verbose: bool
            Print model summary if nonzero.

        Returns
        -------
        model: tensorflow.keras.models.Sequential
            The created sequential ANN model.
        '''
        self._model = tensorflow.keras.models.Sequential()
        for layer_index, num_nodes in enumerate(self._node):
            self._model.add(
                tensorflow.keras.layers.Dense(
                    units=num_nodes,
                    input_dim=self._inputs if layer_index == 0 else None,
                    activation=self._activation[0],
                    kernel_initializer=self._initializer[0],
                    use_bias=self._bias[0]
                )
            )

        if verbose:
            self._model.summary()

        return self._model

    def prepare_input(self, data: np.ndarray) -> np.ndarray:
        '''
        Reshape/prepare the input data array for use with the model.

        Parameters
        ----------
        data: np.ndarray

        Returns
        -------
        out: np.ndarray
            Reshaped input data array (no copy is made only the shape property
            of the data array is adjusted)
        '''
        data.shape = data.shape[:2]
        return data

    def todict(self) -> dict:
        '''
        Export topology to a dict

        Returns
        -------
        data: dict
            A dictionary with all the topology data required to create
            a model.
        '''
        # update configuration
        return {
            'type': self.__class__.__name__,
            'inputs': self._inputs,
            'node': self._node,
            'activation': self._activation,
            'bias': self._bias,
            'initializer': self._initializer,
            'pp_input': self._pp_input.todict(),
            'pp_output': self._pp_output.todict(),
        }

    def save(self, filename: str):
        '''
        Save network topology to a json file.

        Parameters
        ----------
        filename: str or dict
            Name of the target file. Extension is not required.
        '''
        data = self.todict()
        if not filename.endswith(TOPOLOGY_EXTENSION):
            filename = filename + TOPOLOGY_EXTENSION

        with open(filename, 'w') as fid:
            json.dump(data, fid, indent=2)

    @classmethod
    def fromdict(cls, data: dict) -> 'Topology':
        '''
        Load :py:class:`xopto.util.regression.assn.Topology` instance from a
        dict object.

        Parameters
        ----------
        data: dict
            Topology object exported to a dict.

        Returns
        -------
        topology: Topology
            A new instance of :py:class:`xopto.util.regression.ann.Topology`
            created from the data in the dict.
        '''
        data = dict(data)
        type_name = data.pop('type')
        data['pp_input'] = Preprocessor.fromdict(data.pop('pp_input'))
        data['pp_output'] = Preprocessor.fromdict(data.pop('pp_output'))
        return cls(**data)

    @classmethod
    def load(cls, filename: str) -> 'Topology':
        '''
        Load :py:class:`xopto.util.regression.assn.Topology` instance from a
        JSON file.

        Parameters
        ----------
        filename: str
            Source JSON file. The extension .json is added automatically.
        '''
        if not filename.endswith(TOPOLOGY_EXTENSION):
            filename = filename + TOPOLOGY_EXTENSION

        with open(filename, 'r') as fid:
            data = json.load(fid)

        return cls.fromdict(data)

    def __str__(self):
        return 'Topology(inputs={}, node={}, bias={}, activation={}, ' \
                'initializer={}, pp_input={}, pp_output={})'.format(
                    self.inputs, self.node, self.bias, self.activation,
                    self.initializer, self.pp_input, self.pp_output
                )

    def __repr__(self):
        return '{} # {}'.format(self.__str__(), id(self))


class Model:
    def __init__(self, topology: Topology = None, filename: str = None):
        '''
        Creates an ANN model from a topology object or a model file.
        To use an existing model, provide two files with same name and
        extensions .h5 (keras model file) and .json (topology file). To
        create a new model, provide a topology object (an instance of
        Topology).

        Parameters
        ----------
        topology: Topology
            Topology object or None.
        filename: str or (str, str)
            Keras model (.h5) file and topology (.json) file
            (extension not required).
            If list or tuple, the first item should be the model file and
            the second item should be the topology file.
        '''
        self._model_filename = self._topology_filename = None

        if filename is not None:
            if isinstance(filename, (list, tuple)):
                model_filename = str(filename[0])
                topology_filename = str(filename[1])
            else:
                model_filename = topology_filename = str(filename)

            if model_filename.endswith((MODEL_EXTENSION, TOPOLOGY_EXTENSION)):
                model_filename = os.path.splitext(model_filename)[0]
            if topology_filename.endswith((MODEL_EXTENSION, TOPOLOGY_EXTENSION)):
                topology_filename = os.path.splitext(topology_filename)[0]

            model_filename += MODEL_EXTENSION
            topology_filename += TOPOLOGY_EXTENSION

            self._model_filename = model_filename
            self._topology_filename = topology_filename

            model = tensorflow.keras.models.load_model(model_filename)
            with open(topology_filename, 'r') as fid:
                data = json.load(fid)
            topology = Topology.fromdict(data['topology'])

        else:
            model = topology.model()

        self._fit_configuration = {
            'batch_size': [],
            'epochs': [],
            'lr': [],
            'loss': [],
            'save_best': [],
            'validation': [],
            'patience': [],
            'min_delta': [],
            'bestfile': [],
            'shuffle': [],
            'validation_split': [],
        }
        self._fit_results = {
            'loss_history': [],
            'val_loss_history': [],
            'fit_time': []
        }

        self._topology = topology
        self._model = model
        self._train_data = None

    def _get_noutput(self) -> int:
        return self._model.output.shape[-1]
    noutput = property(_get_noutput, None, None, 
                        'The number of outputs.')

    def _get_ninput(self) -> int:
        return self._model.input.shape[-1]
    ninput = property(_get_ninput, None, None, 
                       'The number of inputs.')

    def _get_topology(self) -> Topology:
        return self._topology
    topology = property(_get_topology, None, None,
                        'Topology instance used by this model.')

    def weights(self, layer_index: int) -> List[np.ndarray]:
        '''
        Returns the weights of the selected layer.

        Parameters
        ----------
        layer_index: int
            Zero-based index of the layer.

        Returns
        -------
        weights, bias: List[np.ndarray]
            Weight of the layer connections as a Numpy array.
            Bias terms of the layer as a numpy vector.

        Note
        ----
        Output of a layer is computed as :math:`w \cdot x + b`, where
        :math:`w` is a matrix of weights of the connections, :math:`b` is a
        vector of bias terms and :math:`x` is the input of the layer.
        '''
        return self._model.layers[layer_index].get_weights()

    def activation(self, layer_index: int) -> object:
        '''
        Returns the activation function instance of the layer.

        Parameters
        ----------
        layer_index: int
            Zero-based index of the layer.

        Returns
        -------
        activation: object
            Weight of the layer connections as a Numpy array.
        '''
        return self._model.layers[layer_index].activation

    def _get_layers(self) -> List[tensorflow.keras.layers.Dense]:
        return self._model.layers
    layers = property(_get_layers, None, None, 'List of layers.')

    def summary(self, *args, **kwargs) -> str:
        '''
        Prints a string summary of the model.

        Parameters
        ----------
        args:
            Optional positional arguments passed to the model method.
        kwargs:
            Optional keyword arguments passed to the model method.
        '''
        return self._model.summary(*args, **kwargs)

    def prepare_train(self,
                      train_input: np.ndarray,
                      train_output: np.ndarray,
                      validate_input: np.ndarray = None,
                      validate_output: np.ndarray = None):
        '''
        Training and validation data sets are prepared using the
        pre/post-processors of the input and output variables defined by
        the topology.
        Validation data set is optional and is used to estimate the point
        when the training process can be stopped.

        Parameters
        ----------
        train_input: np.ndarray
            2D training data set of inputs with shape (m, n) where m
            is the number of observations (e.g. wavelengths) and n is the
            number of variables (e.g. source detector separations).
        train_output: np.ndarray
            1D or 2D array of true values with shape (m, k) where m
            is the number of observations (e.g. wavelengths) and k is the
            number of predicted variables (e.g. mua, musr, gamma, delta, ...).
            Usually, k=1 because each parameter is predicted by a dedicated
            network/topology.
        validate_input: np.ndarray
            2D array with input validation data set. The second dimension
            must equal the second dimension of the train_input data set.
        validate_output: np.ndarray
            1D or 2D array of true values of the validation data set. The
            second dimension must equal the second dimension of the
            train_output data set.
        '''
        #####################################################################
        # compute/train the pre/post-processor statistics used for the inputs
        self._topology.pp_input.train(train_input)
        self._topology.pp_output.train(train_output)
        #####################################################################

        # prepare data
        topology = self.topology

        train_input = topology.pp_input.apply(train_input)
        train_output = topology.pp_output.apply(train_output)

        if validate_input is not None:
            validate_input = topology.pp_input.apply(validate_input)
            validate_output = topology.pp_output.apply(validate_output)

        self._train_data = {
            'train_input': train_input,
            'train_output': train_output,
            'validate_input': validate_input,
            'validate_output': validate_output
        }

    def train(self,
              batch_size: int = 32,
              epochs: int =1000,
              optimizer=None,
              lr: float = 1e-3,
              lrcb: Callable[[int, float], float] = None,
              loss: str = 'mse',
              loss_weights=None,
              bestfile: str = None,
              patience: bool = 20, min_delta: float = 0.1e-5,
              shuffle: bool = True,
              validation_split: float = 0.0,
              verbose: bool = False,):
        '''
        Fit/train the ANN model. If a validation data set is passed to the
        prepareTrain method, the training process can be stopped early
        depending on the values of the patience and min_delta parameters.

        Parameters
        ----------
        batch_size: int
            Batch size.
        epochs: int
            Number of epochs.
        verbose: int
            Verbose mode.
        optimizer: str or keras.optimizers.Adam, ...
            One of the Keras optimizers. Default is Adam.
        lr: float
            Default Adam optimizer learning rate. Only used if the optimizer
            parameter is not specified (None).
        lrcb: Callable[[int, float], float]
            A callable that takes takes an epoch index as input
            (integer, indexed from 0) and current learning rate and
            returns a new learning rate as output (float).
        loss: str
            Loss metric.
        loss_weights: List[float] or dict
            Optional list or dictionary specifying scalar coefficients
            (Python floats) to weight the loss contributions of different model
            outputs. If a list, it is expected to have a 1:1 mapping to the
            model's outputs. If a tensor, it is expected to map output names
            (strings) to scalar coefficients.
        bestfile: str
            After each training epoch, the intermediate model (if improvement
            is observed) is optionally saved to the specified file.
        patience:int
            Number of epochs with no improvement after which the training
            will be stopped
        min_delta: float
            Minimum change in the monitored quantity to qualify as an
            improvement, i.e. an absolute change of less than min_delta,
            will not count as an improvement.
        shuffle: bool
            Training data set is randomly shuffled if set to True.
        validation_split: float
            Float between 0 and 1. Fraction of the training data to be used as
            validation data. The model will set apart this fraction of the
            training data, will not train on it, and will evaluate the loss
            and any model metrics on this data at the end of each epoch.
            The validation data is selected from the last samples in the x
            and y data provided, before shuffling. This option is ignored in
            case a dedicated validation set is provided to the prepareTrain
            method.
        validation_freq: int
            Only relevant if validation data is provided. Integer or
            list/tuple/set. If an integer, specifies how many training epochs
            to run before a new validation run is performed,
            e.g. validation_freq=2 runs validation every 2 epochs. If a list,
            tuple, or set, specifies the epochs on which to run validation,
            e.g. validation_freq=[1, 2, 10] runs validation at the end of
            the 1st, 2nd, and 10th epochs.
            Note: Parameter supported only for keras version >= "2.2.4"!
        '''

        validation_split = float(validation_split)
        shuffle = bool(shuffle)

        if shuffle:
            # shuffle the training data set
            n = self._train_data['train_input'].shape[0]
            ind = np.arange(n)
            np.random.shuffle(ind)

            self._train_data['train_input'] = \
                self._train_data['train_input'][ind, ...]
            self._train_data['train_output'] = \
                self._train_data['train_output'][ind, ...]

        for key, value in self._train_data.items():
            if value is not None:
                self._train_data[key] = self._topology.reshapeInput(value)

        validate = False
        if self._train_data['validate_input'] is not None or \
                validation_split > 0.0:
            validate = True

        # create a custom callback monitor
        best_val_loss = float('inf')
        best_val_loss_epoch = 0
        class CustomCallbackMonitor(tensorflow.keras.callbacks.Callback):
            def __init__(self):
                self._best_val_loss = float('inf')
                self._best_val_loss_epoch = 0
                super(CustomCallbackMonitor, self).__init__()

            def on_epoch_end(self, epoch, logs=None):
                if logs is not None:
                    if validate:
                        if logs['val_loss'] < self._best_val_loss:
                            self._best_val_loss = logs['val_loss']
                            self._best_val_loss_epoch = epoch + 1

                        info = '\r  epoch {:d}/{:d}, '\
                               'loss: {:.2e}, val_loss: {:.2e} '\
                               '(best: {:.2e} @epoch {:d})'.format(
                                   epoch + 1, self.params['epochs'],
                                   logs['loss'], logs['val_loss'],
                                   self._best_val_loss,
                                   self._best_val_loss_epoch)
                    else:
                        info = '\r  epoch {:d}/{:d}, loss: {:.2e}'.format(
                            epoch + 1, self.params['epochs'], logs['loss'])
                sys.stdout.write(info)
                sys.stdout.flush()

        custom_monitor = CustomCallbackMonitor()
        callbacks_list = [custom_monitor]

        # learning rate scheduler
        if lrcb is not None:
            lr_callback = tensorflow.keras.callbacks.LearningRateScheduler(lrcb)
            callbacks_list.append(lr_callback)

        print()

        if optimizer is None:
            optimizer =  tensorflow.keras.optimizers.Adam(lr=lr)

        start_training_timestamp = time.perf_counter()
        self._model.compile(optimizer=optimizer, loss=loss,
                            loss_weights=loss_weights)

        if self._train_data['validate_input'] is not None:
            early_stopping = tensorflow.keras.callbacks.EarlyStopping(
                monitor='val_loss',
                min_delta=min_delta,
                patience=patience,
                verbose=False
            )

            callbacks_list.append(early_stopping)
            if bestfile is not None:
                if not bestfile.endswith(MODEL_EXTENSION): 
                    bestfile = bestfile + MODEL_EXTENSION

                checkpoint = tensorflow.keras.callbacks.ModelCheckpoint(
                    bestfile,
                    monitor='val_loss',
                    verbose=False,
                    save_best_only=True,
                    mode='auto'
                )
                callbacks_list.append(checkpoint)

            fit = self._model.fit(
                self._train_data['train_input'],
                self._train_data['train_output'],
                epochs=epochs,
                batch_size=batch_size,
                verbose=verbose,
                validation_data=
                (self._train_data['validate_input'],
                 self._train_data['validate_output']),
                callbacks=callbacks_list,
                validation_split=validation_split,
                shuffle=shuffle,
            )
        else:
            fit = self._model.fit(
                self._train_data['train_input'],
                self._train_data['train_output'],
                epochs=epochs,
                batch_size=batch_size,
                verbose=verbose,
                callbacks=callbacks_list,
                validation_split=validation_split,
                shuffle=shuffle,
            )

        fit_time = time.perf_counter() - start_training_timestamp

        self._fit_configuration['batch_size'].append(batch_size)
        self._fit_configuration['validation_split'].append(validation_split)
        self._fit_configuration['shuffle'].append(shuffle)
        self._fit_configuration['epochs'].append(epochs)
        self._fit_configuration['lr'].append(lr)
        self._fit_configuration['loss'].append(loss)
        self._fit_configuration['bestfile'].append(bestfile)
        if self._train_data['validate_input'] is not None:
            self._fit_configuration['validation'].append('True')
            self._fit_configuration['patience'].append(patience)
            self._fit_configuration['min_delta'].append(min_delta)

        self._fit_results['loss_history'].append(fit.history['loss'])
        if self._train_data['validate_input'] is not None:
            self._fit_results['val_loss_history'].append(
                fit.history['val_loss'])
        self._fit_results['fit_time'].append(fit_time)

    def train_report(self) -> dict:
        '''
        Returns a dict with train statistics and data.

        Returns
        -------
        report: dict
            Train statistics and data.
        '''
        return self._fit_results

    def train_configuration(self) -> dict:
        '''
        Returns a dict with train configuration.

        Returns
        -------
        cfg: dict
            Train configuration.
        '''
        return self._fit_configuration

    def predict(self, input: np.ndarray) -> np.ndarray:
        '''
        Make prediction with the model. Prediction are made in three steps:
        1. Preprocessing (normalization, scaling, ...) of the input data.
        2. Predicting outputs by the model.
        3. Postprocessing (denormalization, descaling, ...) of the predictions.

        Parameters
        ----------
        input: np.ndarray
            2D array of shape (m, n) where m is he number of
            observations (e.g. wavelengths) and n is the number of
            variables (e.g. source-detector separations).

        Returns
        -------
        y: np.ndarray
            Numpy array of predicted values.
        '''
        input = self._topology.prepare_input(input)
        input = self._topology.pp_input.apply(input)
        y = self._model.predict(input)
        y = self._topology.pp_output.undo(y)
        if y.shape[1] == 1:
            y.shape = (y.size,)
        return y

    def benchmark(self, input: np.ndarray, true_value: np.ndarray,
                  threshold: float = None) \
                      -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        '''
        Predictions are evaluated by the RMSE and relative RMSE metrics.
        Optionally, samples for which the true values do not exceed the
        given threshold can be excluded from the relative RMSE estimation.

        Parameters
        ----------
        input: np.ndarray
            Array of input samples stored row-wise.
        true_value: np.ndarray
            Array of true values that will be assessed against the model
            predictions using the RMSE and relative RMSE metrics.
        threshold: float
            Threshold for relative rmse benchmark. Only the true values
            that exceed this threshold are included in the evaluation.

        Returns
        -------
        y: np.ndarray
            Predicted values.
        rmse: float
            Root mean square error
        rrmse: float
            Relative root mean square error
        '''
        predictions = self.predict(input)  # [:, 0]

        rmse = Errors.rmse(true_value, predictions)
        rrmse = Errors.relative_rmse(true_value, predictions, threshold)

        return predictions, rmse, rrmse

    def save(self, filename):
        '''
        Save trained ANN model for later use. The ANN model is saved in the
        Keras native .h5 file format. In addition, a file with the same name
        and .json extension is created. This file holds the corresponding
        topology data (key "topology"), fit configuration
        (key "fit_configuration") and fit results (key "fit_configuration").

        Parameters
        ----------
        filename: str
            Full file name (extension is not required).
        '''
        if filename.endswith((MODEL_EXTENSION, TOPOLOGY_EXTENSION)):
            filename = os.path.splitext(filename)[0]

        self._model_filename = self._topology_filename = filename
        
        model_filename = filename + MODEL_EXTENSION
        topology_filename = filename + TOPOLOGY_EXTENSION
        self._model.save(model_filename)
        data = {
            'topology': self._topology.todict(),
            'fit_configuration': self._fit_configuration,
            'fit_results': self._fit_results,
        }

        with open(topology_filename, 'w') as fid:
            json.dump(data, fid, indent=4)

    @classmethod
    def load(cls, filename: str) -> 'Model':
        '''
        Load model, underlaying topology, fit configuration and fit results
        from an existing '.h5' and '.json' files.

        Returns
        -------
        model: Model
            A :py:class:`xopto.util.regression.ann.Model` instance created
            from the data in the files.
        '''
        return Model(filename=filename)


    def filename(self) -> str:
        '''
        Name of the file including path from which the current model was
        loaded or saved to.

        Returns
        -------
        filename: str
            File name.
        '''
        return self._model_filename


if __name__ == '__main__':
    topology = Topology(node=[40, 20, 10, 5, 1])
    topology.model(verbose=True)
    topology.save('/home/miran/tmp/topology.json')
    t1 = Topology.load('/home/miran/tmp/topology.json')

    d = topology.todict()
    d1 = Topology.fromdict(d)
