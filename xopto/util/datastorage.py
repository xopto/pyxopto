from typing import List, Tuple, Dict
import os
import os.path
import glob
import pickle
import re

import numpy as np


class DatasetFiles:
    @staticmethod
    def analyze_indices(indices: np.ndarray) \
            -> Tuple[List[Tuple[int, int]], List[Tuple[int, int]]]:
        '''
        Analyzes indices as returned by the :py:meth:`Dataset.data_files`
        method.

        Parameters
        ----------
        indices: np.ndarray
            A 2D numpy array of item indices contained within the returned
            files. The indices are organized as
            :py:code:`indices[file_index] => (first, last)`,
            where :py:code:`first` is the index of the first item in the file
            and :py:code:`last` the index of last plus 1 item in the file.

        Returns
        -------
        span: Tuple[int, int]
            Span of data items.
        missing: List[Tuple[int, int]]
            A list of missing data items. Each missing item in the list
            is a tuple (first, n), where first is the index of the first
            missing item and n the number of missing items that follow the
            first missing item.
        overlapping: List[Tuple[int, int]]
            A list of overlapping data items. Each overlapping data item
            in the list is a tuple (first, n), where first is the index
            of the first overlapping item and n the number of overlapping items
            that follow the first overlapping item.
        '''
        missing = []
        overlapping = []
        span = (0, 0)

        if indices.size > 0:
            span = (indices[0, 0], indices[-1, 1])
            delta = indices[1:, 0] - indices[:-1, 1]
            for file_index, item_delta in enumerate(delta):
                if item_delta > 0:
                    start = indices[file_index, 1]
                    missing.append((start, item_delta))
                elif item_delta < 0:
                    last = indices[file_index, 1]
                    overlapping.append((last + item_delta, -item_delta))

        return span, missing, overlapping

    @staticmethod
    def parse_filename(filename: str):
        '''
        Parse a file name and return relevant parts.

        Returns
        -------
        start: int
            Index of the first item in the file.
        stop: int
            Index of the last item + 1 in the file.
        prefix: str
            Data file name prefix.
        '''
        pattern = re.compile('(?P<prefix>.*)_(?P<start>[0-9]+)_(?P<stop>[0-9]+)')

        start = stop = prefix = None

        dir_part, name_part = os.path.split(filename)
        base, ext = os.path.splitext(name_part)
        matching_result = pattern.match(base)
        if matching_result is not None:
            start = int(matching_result['start'])
            stop = int(matching_result['stop'])
            prefix = str(matching_result['prefix'])

        return start, stop, prefix

    def __init__(self, filenames: List[str], indices: np.ndarray = None,
                 prefix: str = None):
        '''
        Creates a dataset filens object.

        Parameters
        ----------
        filenames: List[str]
            A list of dataset files.
        indices: np.ndarray
            A numpy array of item indices stored in each file, one row
            per each file as [start_index, stop_index]. Items with indices
            from and including start_index and up to and including
            stop_index - 1 can be found in the corresponding file.
            If None, the indices are extracted from the filenames.
        prefix: str
            Prefix of the dataset files. Must be the same across all the files.
        '''
        self._filenames = self._prefix = None
        self._span = self._missing = self._overlapping = None

        if indices is None or prefix is None:
            indices = []
            for filename in filenames:
                start, stop, prefix = DatasetFiles.parse_filename(filename)

                if self._prefix is None:
                    self._prefix = prefix
                else:
                    if prefix != self._prefix:
                        raise ValueError('Dataset files must have the same prefix!')

                if start is not None and stop is not None:
                    indices.append([start, stop])
                else:
                    raise ValueError('Failed to parse file name "{}"!'.format(filename))
            indices = np.asarray(indices)
            
        if indices.size > 0:
            sort_indices = np.argsort(indices[:, 0])
            indices[:, 0] = indices[:, 0][sort_indices]
            indices[:, 1] = indices[:, 1][sort_indices]

            filenames = [filenames[index] for index in sort_indices]

        span, missing, overlapping = DatasetFiles.analyse_indices(indices)
        
        self._indices = indices
        self._filenames = filenames
        self._span = span
        self._missing = missing
        self._overlapping = overlapping

    def _get_indices(self) -> np.ndarray:
        return self._indices
    indices = property(
        _get_indices, None, None,
        'A 2D numpy array of item indices contained within the files. '
        'The indices are organized as '
        'indices[file_index] => (first, last),'
        'where first is the index of the first item in '
        'the corresponding file and last the index of the last '
        'plus 1 item in the corresponding file.')

    def _get_prefix(self) -> str:
        return self._prefix
    prefix = property(_get_prefix, None, None,
                      'Prefix of the data file names.')

    def _get_files(self) -> Tuple[str]:
        return tuple(self._filenames)
    files = property(_get_files, None, None,
                     'Ordered list of file names.')

    def _get_span(self) -> List[Tuple[int, int]]:
        return self._span
    span = property(_get_span, None, None,
                    'Index of the first and index '
                    'of the last + 1 item in the files')

    def _get_missing(self) -> List[Tuple[int, int]]:
        return self._missing
    missing = property(
        _get_missing, None, None,
        'A list of missing data items. Each missing item in the list '
        'is a tuple (first, n), where first is the index of the first '
        'missing item and n the number of missing items that follow the '
        'first missing item.')
    
    def _get_overlapping(self) -> List[Tuple[int, int]]:
        return self._overlapping
    overlapping = property(
        _get_overlapping, None, None,
        'A list of overlapping data items. Each overlapping data item '
        'in the list is a tuple (first, n), where first is the index '
        'of the first overlapping item and n the number of overlapping items '
        'that follow the first overlapping item.')

    def __getitem__(self, index) -> str:
        return self._filenames[index]

    def __len__(self) -> int:
        return len(self._filenames)


class Dataset:
    def __init__(self, location: str = None, create: bool = False,
                 overwrite=False):
        '''
        Interface that allows saving and loading datasets to a sequence
        of compressed numpy (.npz) or pickle (.pkl) files.

        Parameters
        ----------
        location: str
            Root directory of the dataset. If location is None, the current
            working directory is used as the dataset location.
        create: bool
            Create location and subdir directories if not exist.
        overwrite: bool
            Overwrite existing files if set to True.
        '''
        if location is None:
            location = os.path.abspath(os.getcwd())

        self._location = str(location)
        self._create = bool(create)
        self._overwrite = bool(overwrite)

        if not os.path.isdir(self._location):
            if not self._create:
                raise ValueError('Dataset location does not exist!')
            else:
                os.makedirs(location)

    def _parse_filenames(self, filenames: List[str]) \
            -> Tuple[List[str], np.ndarray]:
        '''
        Parse and sort the data files in the list.

        Parameters
        ----------
        filenames: List[str]
            A list of data files.

        Returns
        -------
        sorted: List[str]
            Sorted data files.
        indices: np.ndarray
            A 2D numpy array of item indices contained within the returned
            files. The indices are organized as
            :py:code:`indices[file_index] => (first, last)`,
            where :py:code:`first` is the index of the first item in the file
            and :py:code:`last` the index of last plus 1 item in the file.
        '''
        indices = []
        data_files = []
        for index, filename in enumerate(filenames):
            start, stop, _ = DatasetFiles.parse_filename(filename)
            if start is not None and stop is not None:
                indices.append([start, stop])
                data_files.append(filename)

        indices = np.asarray(indices)
        sorted_data_files = []
        if indices.size > 0:
            sort_indices = np.argsort(indices[:, 0])
            indices[:, 0] = indices[:, 0][sort_indices]
            indices[:, 1] = indices[:, 1][sort_indices]

            for index in sort_indices:
                sorted_data_files.append(data_files[index])

        return sorted_data_files, indices

    def files(self, ext: str = None, prefix: str = None,
              subdir: str = None) -> DatasetFiles:
        '''
        Returns an instance of :py:class:`DatasetFiles`.

        Parameters
        ----------
        ext: str
            Data file extension (".npz", ".pkl"). Defaults to ".pkl".
        prefix: str
            File name prefix. Defaults to "data".
        subdir: str
            Subdirectory of the data files.

        Returns
        -------
        files: DatasetFiles
            Object that contains information of the dataset files.
        '''
        if prefix is None:
            prefix = 'data'
        if ext is None:
            ext = '.pkl'

        location = self._location
        if subdir is not None:
            location = os.path.join(location, subdir)

        ext = str(ext).lower()
        if ext not in ('.pkl', '.npz'):
            raise ValueError(
                'Unsupported data file extension"{}"!'.format(ext))

        data_files = glob.glob(
            os.path.join(location, '{:s}_*{:s}'.format(prefix, ext))
        )

        data_files, indices = self._parse_filenames(data_files)

        return DatasetFiles(data_files, indices, prefix)

    
    def find(self, subdir: str = None) -> Tuple[DatasetFiles]:
        '''
        Find data file sequences.
        
        Parameters
        ----------
        subdir: str
            Optional relative path to the data directory.
        ext: str
            Optional file extension (".pkl", ".npz").

        Returns
        -------
        dataset: List[DatasetFiles]
            Found sequences of data files.
        '''
        location = self._location
        if subdir is not None:
            location = os.path.join(location, subdir)

        groups = {}

        for ext in ('.npz', '.pkl'):
            data_files = glob.glob(os.path.join(location, '*' + ext))
            for data_file in data_files:
                start, stop, prefix = DatasetFiles.parse_filename(data_file)
                if start is not None and stop is not None and prefix is not None:
                    group = groups.get(prefix, [])
                    group.append(data_file)
                    groups[prefix] = group
            
        return [DatasetFiles(filenames) for prefix, filenames in groups.items()]

    def load(self, filename: str) -> object:
        '''
        Load data from a .pkl or .npz file.

        Parameters
        ----------
        filename: str
            File from which to load data. File extension will be used
            to resolve the file type.

        Returns
        -------
        data: object
            File content.
        '''
        ext = os.path.splitext(filename)[-1]

        if ext not in ('.pkl', '.npz'):
            raise ValueError(
                'File extension "{}" is not supported!'.format(ext))

        if ext == '.npz':
            data = np.load(filename)
        else:
            with open(filename, 'rb') as fid:
                data = pickle.load(fid)

        return data

    def save(self, data: object, first: int, n: int,
             prefix: str = None, subdir: str = None, ext: str = None,
             silent: bool = False) -> str:
        '''
        Save data items to a compressed structured numpy or pickle file.

        Parameters
        ----------
        data: object
            Data to save.
        first: int
            Index of the first item in the data.
        n: int
            Number of items in the data.
        prefix: str
            File name prefix. Defaults to "data".
        ext: str
            Select file extension/type (".pkl", ".npz"). If saving to NPZ file,
            data that are an instance of dict are expanded. Defaults to ".pkl".
        subdir: str
            Subdirectory for the data files.
        silent: bool
            Does not raise errors if the file exists and cannot be overwritten.

        Returns
        -------
        filename: str
            The name of the saved file.
        '''
        if prefix is None:
            prefix = 'data'
        
        if ext is None:
            ext = '.pkl'

        location = self._location
        if subdir is not None:
            location = os.path.join(location, subdir)

        if not os.path.isdir(location):
            if not self._create:
                raise ValueError('Not allowed to create target '
                                 'directory "{}"!'.format(location))
            else:
                os.makedirs(location, exist_ok=True)

        filename = '{}_{}_{}{}'.format(prefix, first, first + n, ext)
        target = os.path.join(location, filename)
        if os.path.isfile(target) and not self._overwrite:
            if not silent:
                raise ValueError(
                    'Cannot overwrite existing file "{}"!'.format(target))
            return

        if ext == '.npz':
            if isinstance(data, dict):
                np.savez_compressed(target, **data)
            else:
                np.savez_compressed(target, data)
        else:
            with open(target, 'wb') as fid:
                pickle.dump(data, fid)

        return target

    def _get_location(self) -> str:
        return self._location
    location = property(_get_location, None, None, 'Dataset location.')

    def __str__(self):
        return 'Dataset(location="{}", create={}, overwrite={})'.format(
            self._location, self._create, self._overwrite)

    def __repr__(self):
        return '{} # {}'.format(self.__str__(), id(self))

