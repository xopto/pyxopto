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

from typing import Tuple, Callable, Literal

import math
import random
from collections import namedtuple
import json
import traceback

from PySide6.QtWidgets import QApplication, QWidget, QVBoxLayout, QHBoxLayout, \
     QLineEdit, QLabel, QComboBox, QRadioButton, QGroupBox, QGridLayout, \
     QMessageBox, QButtonGroup, QScrollArea, QPushButton, QListWidget, QDialog, \
     QListWidgetItem, QStyle, QFrame, QSizePolicy, \
     QTableWidget, QTableWidgetItem, QHeaderView, \
     QFileDialog
from PySide6.QtCore import Qt, Signal as QSignal, QSize
from PySide6.QtGui import QValidator, QDoubleValidator, QColor, QPixmap, QIcon

from xopto.materials import ri, density
from xopto.materials.ri.base import RefractiveIndex
from xopto.util.suspension import Suspension
from xopto.util.distribution import Normal 

Material = namedtuple('Material', ('id', 'name', 'ri', 'density'))
MATERIAL_POLYSTYRENE = Material(0, 'Polystyrene', ri.polystyrene.default, density.polystyrene.default)
MATERIAL_QUARTZ = Material(1, 'Quartz', ri.glass.fusedsilica.default, 2203.0)
MATERIAL_PMMA = Material(2, 'PMMA', ri.pmma.default, 1185.0)
MATERIAL_CUSTOM = Material(3, 'Custom', None, None)
MATERIAL_WATER = Material(4, 'Water', ri.water.default, density.water.default)

MATERIALS = [
    MATERIAL_POLYSTYRENE, MATERIAL_QUARTZ, MATERIAL_PMMA,
    MATERIAL_CUSTOM, MATERIAL_WATER
]
NAME_TO_MATERIAL = dict([(material.name, material) for material in MATERIALS])

PARTICLE_MATERIALS = [
    MATERIAL_POLYSTYRENE, MATERIAL_QUARTZ, MATERIAL_PMMA, MATERIAL_CUSTOM, MATERIAL_WATER,
    MATERIAL_CUSTOM
]
MEDIUM_MATERIALS = [
    MATERIAL_WATER, MATERIAL_CUSTOM
]

Tolerance = namedtuple('Tolerance', ('id', 'name'))

TOLERANCE_STD = Tolerance(0, 'STD')
TOLERANCE_STD_PERCENT = Tolerance(0, 'STD (%)')

TOLERANCES_DEFAULT = (TOLERANCE_STD, TOLERANCE_STD_PERCENT)


class ErrorMessage(QMessageBox):
    def __init__(self, parent=None):
        super().__init__(parent=parent)
        self.setWindowTitle('Error')
        self.setIcon(QMessageBox.Critical)

    def showMessage(self, error: str, details: str or None = None):
        self.setText(error)
        self.setDetailedText(details if details is not None else '')
        return self.exec()


class WarningConfirmMessage(QMessageBox):
    def __init__(self, parent=None,
                 buttons: Tuple[QMessageBox.StandardButton, ...] or None = None):
        if buttons is None:
            buttons = (self.StandardButton.Yes, self.StandardButton.No)
        super().__init__(parent=parent)
        self.setWindowTitle('Warning')
        self.setIcon(QMessageBox.Warning)
        for button in buttons:
            self.addButton(button)

    def showMessage(self, error: str, details: str or None = None):
        self.setText(error)
        self.setDetailedText(details if details is not None else '')
        return self.exec()


class ValidationError(Exception):
    def __init__(self, msg):
        super().__init__(msg)

    def _get_message(self):
        return self.args[0]

    message = property(_get_message, None ,None)


class Validator:
    @classmethod
    def fromdict(cls, data:dict) -> 'Validator':
        data = dict(data)
        t = data.pop('type')
        if isinstance(t, str):
            t = globals().get(t)
        if not issubclass(t, cls):
            raise TypeError('Expected data of instance {} but got {}!'.format(
                cls.__name__, t.__name__))

        return t(**data)


class RangeValidator(Validator):
    def __init__(self, low: float, high: float, finite: bool = True):
        self.low = float(low)
        self.high = float(high)
        self.finite = bool(finite)

    def todict(self) -> dict:
        return {
            'type': self.__class__.__name__,
            'low': self.low, 'high': self.high, 'finite': self.finite
        }

    def __call__(self, value: float) -> bool:

        if self.finite and not math.isfinite(value):
            raise ValidationError('Value must be finite!')

        if not (self.low <= value <= self.high):
            raise ValidationError(
                'Value must be in range from {} to {}!'.format
                (self.low, self.high))

    def __repr__(self):
        return 'RangeValidator(low={}, high={}, finite={})'.format(
            self.low, self.high, self.finite)

    def __str__(self):
        return '{} # {}'.format(self.__repr__(), id(self))


class NonzeroPositiveValidator(Validator):
    def __init__(self, finite: bool = True):
        self.finite = bool(finite)

    def __call__(self, value: float):
        if self.finite and not math.isfinite(value):
            raise ValidationError('Value must be finite!')

        if value <= 0.0:
            raise ValidationError('Value must be positive and nonzero!')

    def todict(self) -> dict:
        return {'type': self.__class__.__name__, 'finite': self.finite}

    def __repr__(self):
        return 'NonzeroPositiveValidator(finite={})'.format(self.finite)

    def __str__(self):
        return '{} # {}'.format(self.__repr__(), id(self))

class NonnegativeValidator(Validator):
    def __init__(self, finite: bool = True):
        self.finite = bool(finite)

    def __call__(self, value: float):
        if self.finite and not math.isfinite(value):
            raise ValidationError('Value must be finite!')

        if value < 0.0:
            raise ValidationError('Value must be nonnegative!')

    def todict(self) -> dict:
        return {'type': self.__class__.__name__, 'finite': self.finite}

    def __repr__(self):
        return 'NonzeroPositiveValidator(finite={})'.format(self.finite)

    def __str__(self):
        return '{} # {}'.format(self.__repr__(), id(self))


class NormalParameter:
    def __init__(self, value, std: float = 0.0,
                 scifactor: float = 1.0, units: str or None = None):

        self._value = float(value)
        self._std = float(std)
        self._units = units
        self._scifactor = float(scifactor)

    def _get_sci_factor(self) -> float:
        return self._scifactor
    def _set_sci_factor(self, factor: float):
        self._scifactor = float(factor)
    sciFactor = property(_get_sci_factor, _set_sci_factor, None,
                         'SCI => value factor')

    def _get_value(self) -> float:
        return self._value
    def _set_value(self, value: float):
        self._value = value
    value = property(_get_value, _set_value, None, 'Value.')

    def _get_sci_value(self) -> float:
        return self._value/self._scifactor
    def _set_sci_value(self, value: float):
        self._value = value*self._scifactor
    sciValue = property(_get_sci_value, _set_sci_value, None, 'Sci value.')

    def _get_units(self) -> str or None:
        return self._units
    def _set_units(self, units: str or None):
        self._units = units
    units = property(_get_units, _set_units, None, 'Value unists.')

    def _get_std(self) -> float:
        return self._std
    std = property(_get_std, None, None, 'Value SCI STD.')

    def _get_sci_std(self) -> float:
        return self._std/self._scifactor
    def _set_sci_std(self, sci: float):
        self._std = sci/self._scifactor
    sciStd = property(_get_sci_std, _set_sci_std, None, 'Value SCI STD.')

    def sample(self, clip: float = 5.0) -> float:
        while abs(self._value - sample) <= clip*self._std:
            sample = random.normalvariate(self._value, self._std)
        return sample

    def sampleSci(self, clip: float = 5.0) -> float:
        return self.sample(clip=clip)/self.sciFactor

    def todict(self) -> dict:
        return {
            'type': self.__class__.__name__,
            'value': self.value, 'std': self.std,
            'units': self.units,
            'scifactor': self.sciFactor
        }

    def tojson(self) -> str:
        return json.dumps(self.todict())

    @classmethod
    def fromdict(cls, data: dict) -> 'NormalParameter':
        data = dict(data)
        t = data.pop('type')
        if isinstance(t, str):
            if t != cls.__name__:
                raise TypeError(
                    'Expected data for type {} but got {}!'.format(
                        cls.__name__, t))
            t = cls
        if t != cls:
            raise TypeError(
                'Expected data of type {} but got {}!'.format(cls, t))

        return cls(**data)

    def update(self, data: dict):
        self.sciFactor = data.get('scifactor', self.sciFactor)
        self.value = data.get('value', self.data)
        self.std = data.get('std', self.std)
        self.units = data.get('units', self.units)
        
    def __str__(self):
        return 'NormalParameter(value={}, std={}, scifactor={}, ' \
                'units="{}")'.format(self.value, self.std,
                                     self.scifactor, self.units)


class NormalParameterWidget(QWidget):
    TOLERANCE_STD = 0
    TOLERANCE_STD_PERCENT = 1

    def __init__(self, name='Parameter', units='', sci_factor=1.0,
                 parent = None):
        '''
        Parameters
        ----------
        name: str
            Name of the parameter.
        units: str
            Units of the parametr.
        sci_factor:
            Factor that converts the native value of parameter to a value
            in SCI units.

        '''
        super().__init__(parent)

        self._tolerances = TOLERANCES_DEFAULT

        grid_layout = QGridLayout()

        self.errorDialog = ErrorMessage(self)
        self.errorDialog.setModal(True)

        self.parameterGroupBox = QGroupBox()
        self.parameterGroupBox.setTitle(name)

        # units
        self.unitsLabel = QLabel()
        if units:
            self.unitsLabel.setText(units)
        self.unitsSciFactor = float(sci_factor)

        # Reference value
        self.valueLabel = QLabel()
        self.valueLabel.setText('Value')
        self.valueLineEdit = QLineEdit()
        self.valueValidator = None

        # Reference value tolerance
        self.valueToleranceLabel = QLabel()
        self.valueToleranceLabel.setText('Tolerance')
        self.valueToleranceLineEdit = QLineEdit()
        self.valueToleranceLineEdit.setText('0.0')
        self.valueToleranceLabel.setBuddy(self.valueToleranceLineEdit)
        self.valueToleranceValidator = None

        # Tolerance kind ... STD or STD %
        self.valueToleranceKindLable = QLabel()
        self.valueToleranceKindLable.setText('Tolerance kind')
        self.valueToleranceKindComboBox = QComboBox()
        self.valueToleranceKindLable.setBuddy(self.valueToleranceKindComboBox)
        self.valueToleranceKindComboBox.addItems(
            [t.name for t in self._tolerances])
        self.valueToleranceKindComboBox.currentIndexChanged.connect(
            lambda ind: self.setValueToleranceKind(self._tolerances[ind]))
        self.valueToleranceKindComboBox.setCurrentIndex(0)
        self.valueLabel.setBuddy(self.valueToleranceKindComboBox)

        grid_layout.addWidget(self.valueLabel, 0, 0)
        grid_layout.addWidget(self.valueToleranceLabel, 0, 1)
        grid_layout.addWidget(self.valueToleranceKindLable, 0, 2)

        grid_layout.addWidget(self.valueLineEdit, 1, 0)
        grid_layout.addWidget(self.valueToleranceLineEdit, 1, 1)
        grid_layout.addWidget(self.valueToleranceKindComboBox, 1, 2)

        grid_layout.addWidget(self.unitsLabel, 1, 3)

        self.parameterGroupBox.setLayout(grid_layout)
        main_layout = QVBoxLayout()
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.addWidget(self.parameterGroupBox)

        self.setLayout(main_layout)

    def todict(self) -> dict:
        valueValidators = valueToleranceValidators = []
        if self.valueValidator is not None:
            valueValidators = [v.todict() for v in self.valueValidator]
        if self.valueToleranceValidator is not None:
            valueToleranceValidators = \
                [v.todict() for v in self.valueToleranceValidator]
        return {
            'type': self.__class__.__name__,
            'name': self.name(),
            'units': self.units(),
            'sciFactor': self.sciFactor(),
            'value': self.value(),
            'valueValidators': valueValidators,
            'tolerance': self.valueTolerance(),
            'toleranceVisible': self.toleranceVisible(),
            'toleranceKind': self.valueToleranceKind(),
            'valueToleranceValidators': valueToleranceValidators
        }

    def fromdict(self, data: dict):
        t = data.get('type')
        if isinstance(t, str):
            if t != self.__class__.__name__:
                raise TypeError(
                    'Expected data for type {} but got {}'.format(
                        self.__class__.__name__, t))

        self.setName(data.get('name'))
        self.setUnits(data.get('units'))
        self.setSciFactor(data.get('sciFactor'))
        self.setValue(data.get('value'))
        self.setValueValidator([Validator.fromdict(validator)
                                for validator in data.get('valueValidators')])

        self.setValueTolerance(data.get('tolerance'))
        self.setValueToleranceKind(Tolerance(*data.get('toleranceKind')))
        self.setValueToleranceValidator(
            [Validator.fromdict(validator)
             for validator in data.get('valueToleranceValidators')])
        self.setTolerancesVisible(data.get('toleranceVisible'))

    def setValueValidator(self, validator: Validator or Tuple[Validator, ...]):
        '''
        Set one or multiple validators for the value of the parameter.
        Validators are applied to value in native units.

        Parameters
        ----------
        validator:  Validator or Tuple[Validator, ...]
            Validator(s) for the parameter value in native units.
        '''
        if isinstance(validator, Validator):
            validator = (validator,)
        self.valueValidator = validator

    def setValueToleranceValidator(
            self, validator: Validator or Tuple[Validator, ...]):
        '''
        Set one or multiple validators for the tolerance of the parameter.
        Validators are applied to the tolerance value in native units.

        Parameters
        ----------
        validator:  Validator or Tuple[Validator, ...]
            Validator(s) for the parameter tolerance in native units.
        '''
        if isinstance(validator, Validator):
            validator = (validator,)
        self.valueToleranceValidator = validator

    def focusValue(self):
        '''
        Move the cursor and focus to the value widget.
        '''
        self.activateWindow()
        self.valueLineEdit.setFocus()

    def focusValueTolerance(self):
        '''
        Move the cursor and focus to the value tolerance widget.
        '''
        self.activateWindow()
        self.valueToleranceLineEdit.setFocus()

    def validateValue(self, prefix: str = None, silent: bool = False) -> bool:
        '''
        Validate the parameter value. Return True if value is defined and valid.

        Parameters
        ----------
        prefix: str
            Prefix included in the error messages.
        silent: bool
            If True, no error dialogs are shown.

        Returns
        -------
        valid: bool
            True if the value is defined and valid.
        '''
        if prefix is not None:
            prefix = prefix.strip(' :') + ':'
        else:
            prefix = ''

        text_value = self.valueLineEdit.text()
        if not text_value:
            if not silent:
                self.errorDialog.showMessage(
                    '{}{}: Enter a value!'.format(prefix, self.name()))
                self.activateWindow()
                self.valueLineEdit.setFocus()
            return False
        value = None
        try:
            value = float(text_value)
        except (ValueError, TypeError):
            pass
        if value is None:
            if not silent:
                self.errorDialog.showMessage(
                    '{}{}: Enter a valid numeric value!'.format(
                        prefix, self.name()))
                self.activateWindow()
                self.valueLineEdit.setFocus()
            return False

        if self.valueValidator:
            for validator in self.valueValidator:
                message = None
                try:
                    validator(value)
                except ValidationError as err:
                    message = err.message

                if message is not None:
                    if not silent:
                        self.errorDialog.showMessage(
                            '{}{}: {}'.format(prefix, self.name(), message))
                        self.activateWindow()
                        self.valueLineEdit.setFocus()
                    return False
                
        return True

    def validateValueTolerance(self, prefix: str = None,
                               silent: bool = False) -> bool:
        '''
        Validate the value tolerance and return True if a tolerance is defined
        and is valid.

        Parameters
        ----------
        prefix: str
            Prefix included in the error messages.
        silent: bool
            If True, no error dialogs are shown.

        Returns
        -------
        valid: bool
            True if the tolerance value is defined and valid.
        '''
        if prefix is not None:
            prefix = prefix.strip(' :') + ':'
        else:
            prefix = ''

        text_value_tol = self.valueToleranceLineEdit.text()
        if not text_value_tol:
            if not silent:
                self.errorDialog.showMessage(
                    '{}{}: Enter a tolerance!'.format(prefix, self.name()))
                self.activateWindow()
                self.valueToleranceLineEdit.setFocus()
            return False
        value = None
        try:
            value = float(text_value_tol)
        except (ValueError, TypeError):
            pass
        if value is None:
            if not silent:
                self.errorDialog.showMessage(
                    '{}{}: Enter a valid numeric tolerance!'.format(
                        prefix, self.name()))
                self.activateWindow()
                self.valueToleranceLineEdit.setFocus()
            return False
        if value < 0:
            if not silent:
                self.errorDialog.showMessage(
                    '{}{}: Tolerance must be a non-negative value!'.format(
                        prefix, self.name()))
                self.activateWindow()
                self.valueToleranceLineEdit.setFocus()
            return False

        if self.valueToleranceValidator:
            for validator in self.valueToleranceValidator:
                message = None
                try:
                    validator(value)
                except ValidationError as err:
                    message = err.message

                if message is not None:
                    if not silent:
                        self.errorDialog.showMessage('{}{}: {}'.format(prefix, self.name(), message))
                        self.activateWindow()
                        self.valueToleranceLineEdit.setFocus()
                    return False
    
        return True

    def validate(self, prefix: str = None, silent: bool = False) -> bool:
        '''
        Validate the value and tolerance.

        Parameters
        ----------
        prefix: str
            Prefix included in the error messages.
        silent: bool
            If True, no error dialogs are showm.

        Returns
        -------
        valid: bool
            True if the value and its tolerance are defined and valid.
        '''
        if not self.validateValue(prefix, silent=silent):
            return False
        if not self.validateValueTolerance(prefix, silent=silent):
            return False

        return True

    def setName(self, name: str):
        '''
        Set the name of the parameter.

        Parameters
        ----------
        name: str
            New name for the parameter.
        '''
        self.parameterGroupBox.setTitle(name)

    def name(self) -> str:
        '''
        Returns
        -------
        name: str
            Current name of the parameter
        '''
        return self.parameterGroupBox.title()

    def setUnits(self, label: str):
        '''
        Set the units of the parameter and the related coversion factor
        from native value to corresponding SCI value.

        Parameters
        ----------
        label: str
            Label of the unit.
        '''
        self.unitsLabel.setText(label)

    def units(self) -> str:
        '''
        Returns
        -------
        label: str
            Label of the unit.
        '''
        return self.unitsLabel.text()

    def sciFactor(self) -> float:
        '''
        Returns
        -------
        factor: float
            Factor that converts a SCI units ta native value.
        '''
        return self.unitsSciFactor
    
    def setSciFactor(self, factor: float):
        '''
        Set the SCI factor of the unit.

        Parameters
        ----------
        factor: float
            Converts SCI value to native value.
        '''
        self.unitsSciFactor = float(factor)

    def valueToleranceKind(self) -> Tolerance or None:
        '''
        Returns the currently selected value tolerance kind.

        Returns
        -------
        tol: Tolerance
            Tolerance kind.
        '''
        ind = self.valueToleranceKindComboBox.currentIndex()
        if ind is not None:
            return self._tolerances[ind]

    def setValueToleranceKind(self, kind: Tolerance):
        '''
        Set the tolerance kind to the specified value.

        Parameters
        ----------
        kind: Tolerance
            Target tolerance kind.
        '''
        if kind not in self._tolerances:
            raise ValidationError('Tolerance kind {} not supported!'.format(kind))
        self.valueToleranceKindComboBox.setCurrentIndex(self._tolerances.index(kind))

    def toleranceVisible(self) -> bool:
        '''
        Returns
        -------
        visible: bool
            True if tolerances are visible, else False.
        '''
        return self.valueToleranceLineEdit.isVisibleTo(self)

    def setTolerancesVisible(self, visible: bool):
        '''
        Sets the visible state of the tolerances.

        Parameters
        ----------
        visible: bool
            Visible state of the tolerances.
        '''
        self.valueToleranceLineEdit.setVisible(visible)
        self.valueToleranceLabel.setVisible(visible)
        self.valueToleranceKindComboBox.setVisible(visible)
        self.valueToleranceKindLable.setVisible(visible)

    def value(self) -> float:
        '''
        Returns
        -------
        value: float
            Value of the parameter in native units or None if the value
            is undefined.
        '''
        if not self.validateValue(silent=True):
            return None

        value = None
        try:
            value = float(self.valueLineEdit.text())
        except (ValueError, TypeError):
            pass

        return value

    def sciValue(self) -> float:
        '''
        Returns
        -------
        value: float
            Parameter value if SCI units or None if the value is undefined.
        '''
        value = self.value()
        if value is not None:
            value = value*self.unitsSciFactor

        return value

    def setValue(self, value: float or None):
        '''
        Set the parameter value in native units.

        Parameters
        ----------
        value: float
            Target value of the parameter in native units.
        '''
        if value is not None:
            self.valueLineEdit.setText('{:f}'.format(value))
        else:
            self.valueLineEdit.setText('')

    def setSciValue(self, value: float or None):
        '''
        Set the parameter value in SCI units.

        Parameters
        ----------
        value: float
            Target value of the parameter in SCI units.
        '''
        if value is not None:
            value = value/self.unitsSciFactor
        self.setValue(value)

    def valueTolerance(self) -> float or None:
        '''
        Get the tolerance of the parameter in native units

        Returns
        -------
        tol: float
            Tolerance of the parameter value in native units
        '''
        if not self.validateValueTolerance(silent=True):
            return None

        if self.valueToleranceKind() == TOLERANCE_STD_PERCENT:
            if not self.validateValue(silent=True):
                return False

        tolerance = None
        try:
            tolerance = float(self.valueToleranceLineEdit.text())
            if self.valueToleranceKind() == TOLERANCE_STD_PERCENT:
                tolerance = tolerance*self.value()/100.0
        except (ValueError, TypeError):
            pass

        return tolerance

    def setValueTolerance(self, value: float):
        '''
        Set the tolerance of the parameter in native units

        Parameters
        ----------
        tol: float
            Tolerance of the parameter value in native units
        '''
        text = ''
        if value is not None:
            if self.valueToleranceKind() == TOLERANCE_STD:
                text = str(value)
            else:
                text = str(100.0*value/self.value())
        self.valueToleranceLineEdit.setText(text)

    def setSciValueTolerance(self, value: float or None):
        '''
        Set the tolerance of the parameter in SCI units

        Parameters
        ----------
        tol: float
            Tolerance of the parameter value in SCI units
        '''
        if value is not None:
            value = value*self.sciFactor()
        self.setValueTolerance(value)

    def sciValueTolerance(self) -> float or None:
        '''
        Get the tolerance of the parameter in SCI units

        Returns
        -------
        tol: float
            Tolerance of the parameter value in SCI units or None if a
            value is not defined.
        '''
        tolerance = None
        try:
            tolerance = float(self.valueToleranceLineEdit.text())
            tolerance = tolerance*self.unitsSciFactor
            if self.valueToleranceKind() == TOLERANCE_STD_PERCENT:
                tolerance = tolerance*self.sciValue()/100.0
        except (ValueError, TypeError):
            pass

        return tolerance

    def parameter(self) -> NormalParameter:
        '''
        Create and return an instance of NormalParameter.
        '''
        return NormalParameter(value=self.value(), std=self.valueTolerance(),
                               scifactor=self.sciFactor(), units=self.units())


class ParticleSizeWidget(QWidget):
    def __init__(self, title: str or None = None, parent=None):
        super().__init__(parent)

        if title is None:
            title = 'Particle size'

        self.groupBox = QGroupBox()
        self.groupBox.setTitle(title)

        self.diameterParameter = NormalParameterWidget('Diameter', 'μm', 1e-6)
        self.diameterParameter.layout().setContentsMargins(0, 0, 0, 0)
        self.diameterParameter.setValueValidator(
            (NonzeroPositiveValidator(), RangeValidator(0.0, 100.0))
        )

        self.diameterStdParameter = NormalParameterWidget('Diameter STD', 'μm', 1e-6)
        self.diameterStdParameter.layout().setContentsMargins(0, 0, 0, 0)

        layout = QVBoxLayout()
        layout.addWidget(self.diameterParameter)
        layout.addWidget(self.diameterStdParameter)
        self.groupBox.setLayout(layout)

        main_layout = QVBoxLayout()
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.addWidget(self.groupBox)

        self.setLayout(main_layout)

    def title(self) -> str:
        return self.groupBox.title()

    def setTitle(self, title: str):
        self.groupBox.setTitle(title)

    def validate(self, silent: bool = False) -> bool:
        if not self.diameterParameter.validate(self.groupBox.title(), silent=silent):
            return False
        if not self.diameterStdParameter.validate(self.groupBox.title(), silent=silent):
            return False

        d, d_std = self.sciDiameter(), self.sciDiameterStd()
        if d - 5*d_std <= 0.0:
            if not silent:
                self.diameterStdParameter.errorDialog.showMessage(
                    '{}:{}: Diameter STD is too large and yields negative '
                    'diameter values with 5 STD clipping!'.format(
                        self.groupBox.title(),
                        self.diameterStdParameter.name()))
                self.diameterStdParameter.focusValue()
            return False

        return True

    def diameter(self) -> float:
        return self.diameterParameter.value()

    def diameterTolerance(self) -> float:
        return self.diameterParameter.valueTolerance()

    def sciDiameter(self) -> float:
        return self.diameterParameter.sciValue()

    def sciDiameterTolerance(self) -> float:
        return self.diameterParameter.sciValueTolerance()

    def diameterStd(self) -> float:
        return self.diameterStdParameter.value()

    def diameterStdTolerance(self) -> float:
        return self.diameterStdParameter.valueTolerance()

    def sciDiameterStd(self) -> float:
        return self.diameterStdParameter.sciValue()

    def sciDiameterStdTolerance(self) -> float:
        return self.diameterStdParameter.sciValueTolerance()

    def todict(self) -> dict:
        return {
            'type': self.__class__.__name__,
            'title': self.title(),
            'diameter': self.diameterParameter.todict(),
            'diameterStd': self.diameterStdParameter.todict(),
        }

    def fromdict(self, data: dict):
        t = data.get('type')
        if t != self.__class__.__name__:
            raise TypeError(
                'Expected data of type {} but got {}!'.format(
                    self.__class__.__name__, t))
        self.setTitle(data.get('title'))
        self.diameterParameter.fromdict(data.get('diameter'))
        self.diameterStdParameter.fromdict(data.get('diameterStd'))


class DensityWidget(QWidget):
    def __init__(self, materials: Tuple[Material, ...],
                 title: str or None = None, parent=None):
        super().__init__(parent)

        if title is None:
            title = 'Particle size'

        self.groupBox = QGroupBox()
        self.groupBox.setTitle(title)
        self._materials = materials

        self.materialComboBox = QComboBox()
        self.materialComboBox.addItems([m.name for m in materials])
        self.materialComboBox.setCurrentIndex(0)
        self.materialComboBox.currentIndexChanged.connect(
            lambda index: self.setMaterial(materials[index]))

        self.customDensityValue = NormalParameterWidget('Custom value', 'kg/m<sup>3</sup>', 1.0)
        self.customDensityValue.setValueValidator(NonzeroPositiveValidator())
        self.customDensityValue.layout().setContentsMargins(0, 0, 0, 0)
        self.customDensityValue.setTolerancesVisible(False)

        self.setMaterial(materials[0])

        layout = QVBoxLayout()
        layout.addWidget(self.materialComboBox)
        layout.addWidget(self.customDensityValue)

        self.groupBox.setLayout(layout)

        main_layout = QVBoxLayout()
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.addWidget(self.groupBox)
        self.setLayout(main_layout)

    def setTitle(self, title: str):
        self.groupBox.setTitle(title)

    def title(self) -> str:
        return self.groupBox.title()

    def material(self) -> Material:
        ind = self.materialComboBox.currentIndex()
        if ind is not None:
            return self._materials[ind]

    def setMaterial(self, material: Material or str):
        if isinstance(material, str):
            material = NAME_TO_MATERIAL[material]
        self.materialComboBox.setCurrentIndex(self._materials.index(material))
        self.customDensityValue.setEnabled(material == MATERIAL_CUSTOM)

    def validate(self, silent: bool = False) -> bool:
        if self.material() == self._materials:
            return self.customDensityValue.validate(
                self.groupBox.title(), silent=silent)
        return True

    def density(self) -> float or density.Density or None:
        material = self.material()
        if material == MATERIAL_CUSTOM:
            return self.customDensityValue.value()
        else:
            return material.density

    def sciDensity(self) -> float or density.Density or None:
        material = self.material()
        if material == MATERIAL_CUSTOM:
            return self.customDensityValue.sciValue()
        else:
            return material.density

    def todict(self) -> dict:
        return {
            'type': self.__class__.__name__,
            'title': self.title(),
            'material': self.material().name,
            'materials': [material.name for material in self._materials],
            'customValue': self.customDensityValue.todict(),
        }


    def fromdict(self, data: dict):
        t = data.get('type')
        if t != self.__class__.__name__:
            raise TypeError(
                'Expected data of type {} but got {}!'.format(
                    self.__class__.__name__, t))
        self.setTitle(data.get('title'))
        self.setMaterial(NAME_TO_MATERIAL[data.get('material')])
        self.customDensityValue.fromdict(data.get('customValue'))


class RefractiveIndexWidget(QWidget):
    def __init__(self, materials = Tuple[Material, ...],
                 title: str or None = None, parent=None):
        super().__init__(parent)

        if title is None:
            title = 'Refractive index'

        self._materials = materials

        self.groupBox = QGroupBox()
        self.groupBox.setTitle(title)

        self.materialComboBox = QComboBox()
        self.materialComboBox.addItems([m.name for m in materials])
        self.materialComboBox.setCurrentIndex(0)
        self.materialComboBox.currentIndexChanged.connect(
            lambda index: self.setMaterial(materials[index]))

        self.customRiValue = NormalParameterWidget('Custom value')
        self.customRiValue.setValueValidator(RangeValidator(1.0, float('inf')))
        self.customRiValue.layout().setContentsMargins(0, 0, 0, 0)
        self.customRiValue.setTolerancesVisible(False)

        self.setMaterial(materials[0])

        layout = QVBoxLayout()
        layout.addWidget(self.materialComboBox)
        layout.addWidget(self.customRiValue)

        self.groupBox.setLayout(layout)

        main_layout = QVBoxLayout()
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.addWidget(self.groupBox)
        self.setLayout(main_layout)

    def setTitle(self, title: str):
        self.groupBox.setTitle(title)

    def title(self) -> str:
        return self.groupBox.title()

    def material(self) -> Material:
        ind = self.materialComboBox.currentIndex()
        if ind is not None:
            return self._materials[ind]

    def setMaterial(self, material: Material or str):
        if isinstance(material, str):
            material = NAME_TO_MATERIAL[material]
        self.materialComboBox.setCurrentIndex(self._materials.index(material))
        self.customRiValue.setEnabled(material == MATERIAL_CUSTOM)

    def validate(self, silent: bool = False) -> bool:
        if self.material() == MATERIAL_CUSTOM:
            return self.customRiValue.validate(
                self.groupBox.title(), silent=silent)
        return True

    def refractiveIndex(self) -> float or RefractiveIndex:
        material = self.material()
        if material == MATERIAL_CUSTOM:
            return self.customRiValue.value()
        else:
            return material.ri

    def sciRefractiveIndex(self) -> float or RefractiveIndex:
        material = self.material()
        if material == MATERIAL_CUSTOM:
            return self.customRiValue.sciValue()
        else:
            return material.ri

    def todict(self) -> dict:
        return {
            'type': self.__class__.__name__,
            'title': self.title(),
            'material': self.material().name,
            'materials': [material.name for material in self._materials],
            'customValue': self.customRiValue.todict(),
        }

    def fromdict(self, data: dict):
        t = data.get('type')
        if t != self.__class__.__name__:
            raise TypeError(
                'Expected data of type {} but got {}!'.format(
                    self.__class__.__name__, t))
        self.setTitle(data.get('title'))
        self.setMaterial(NAME_TO_MATERIAL[data.get('material')])
        self.customRiValue.fromdict(data.get('customValue'))


class ScatteringWidget(QWidget):
    SOLID_CONTENT = 0
    MUSR_WAVELENGTH = 1
    MUS_WAVELENGTH = 2

    def __init__(self, title: str = None, parent=None):
        super().__init__(parent)

        if title is None:
            title = 'Scattering properties'

        self.errorDialog = ErrorMessage(self)
        self.errorDialog.setModal(True)

        self.groupBox = QGroupBox()
        self.groupBox.setTitle(title)

        self.definitionTypeLabel = QLabel()
        self.definitionTypeLabel.setText('Value type')
        self.definitionTypeComboBox = QComboBox()
        self.definitionTypeComboBox.addItems(('Solid content', 'Reduced scattering coefficient', 'Scattering coefficient'))
        self.definitionTypeComboBox.currentIndexChanged.connect(self.setValueType)
        self.definitionTypeLabel.setBuddy(self.definitionTypeComboBox)

        self.solidContentParameter = NormalParameterWidget('Solid content', 'wt %', 1.0)
        self.solidContentParameter.setValueValidator(NonzeroPositiveValidator())
        self.musParameter = NormalParameterWidget('Scattering coefficient', '1/cm', 1e2)
        self.musParameter.setValueValidator(NonzeroPositiveValidator())
        self.musrParameter = NormalParameterWidget('Reduced scattering coefficient', '1/cm', 1e2)
        self.musrParameter.setValueValidator(NonzeroPositiveValidator())
        self.wavelengthParameter = NormalParameterWidget('Wavelength', 'nm', 1e-9)
        self.wavelengthParameter.setTolerancesVisible(False)
        self.wavelengthParameter.setValueValidator(NonzeroPositiveValidator())

        self.solidContentParameter.layout().setContentsMargins(0, 0, 0, 0)
        self.musrParameter.layout().setContentsMargins(0, 0, 0, 0)
        self.musParameter.layout().setContentsMargins(0, 0, 0, 0)
        self.wavelengthParameter.layout().setContentsMargins(0, 0, 0, 0)

        self.setValueType(self.SOLID_CONTENT)

        layout = QVBoxLayout()
        layout.addWidget(self.definitionTypeLabel)
        layout.addWidget(self.definitionTypeComboBox)
        layout.addWidget(self.solidContentParameter)
        layout.addWidget(self.musrParameter)
        layout.addWidget(self.musParameter)
        layout.addWidget(self.wavelengthParameter)

        self.groupBox.setLayout(layout)

        main_layout = QVBoxLayout()
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.addWidget(self.groupBox)
        self.setLayout(main_layout)

    def setTitle(self, title: str):
        self.groupBox.setTitle(title)

    def title(self) -> str:
        return self.groupBox.title()

    def validate(self, silent: bool = False):
        vt = self.valueType()
        if vt == self.SOLID_CONTENT:
            return self.validateSolidContent(self.groupBox.title(), silent=silent)
        elif vt == self.MUS_WAVELENGTH:
            if not self.validateMus(self.groupBox.title(), silent=silent):
                return False
            if not self.validateWavelength(self.groupBox.title(), silent=silent):
                return False
        elif vt == self.MUSR_WAVELENGTH:
            if not self.validateMusr(self.groupBox.title(), silent=silent):
                return False
            if not self.validateWavelength(self.groupBox.title(), silent=silent):
                return False

        return True

    def validateSolidContent(self, prefix: str = None, silent: bool = False) -> bool:
        if prefix is None:
            prefix = self.groupBox.title()
        return self.solidContentParameter.validate(prefix, silent=silent)

    def solidContent(self) -> float:
        return self.solidContentParameter.value()

    def sciSolidContent(self) -> float:
        return self.solidContentParameter.sciValue()

    def solidContentTolerance(self) -> float:
        return self.solidContentParameter.valueTolerance()

    def sciSolidContentTolerance(self) -> float:
        return self.solidContentParameter.sciValueTolerance()

    def validateMus(self, prefix: str = None, silent: bool = False) -> bool:
        if prefix is None:
            prefix = self.groupBox.title()
        return self.musParameter.validate(prefix, silent=silent)

    def mus(self) -> float:
        return self.musParameter.value()

    def sciMus(self) -> float:
        return self.musParameter.sciValue()

    def musTolerance(self) -> float:
        return self.musParameter.valueTolerance()

    def sciMusTolerance(self) -> float:
        return self.musParameter.sciValueTolerance()

    def validateMusr(self, prefix: str = None, silent: bool = False) -> bool:
        if prefix is None:
            prefix = self.groupBox.title()
        return self.musrParameter.validate(prefix, silent=silent)

    def musr(self) -> float:
        return self.musrParameter.value()

    def sciMusr(self) -> float:
        return self.musrParameter.sciValue()

    def musrTolerance(self) -> float:
        return self.musrParameter.valueTolerance()

    def sciMusrTolerance(self) -> float:
        return self.musrParameter.sciValueTolerance()

    def validateWavelength(self, prefix: str = None, silent: bool = False) -> bool:
        if prefix is None:
            prefix = self.groupBox.title()
        return self.wavelengthParameter.validate(prefix, silent=silent)

    def wavelength(self) -> float:
        return self.wavelengthParameter.value()

    def sciWavelength(self) -> float:
        return self.wavelengthParameter.sciValue()

    def setValueType(self, t: int):
        enabled =(True, False, False, False)
        if t == self.SOLID_CONTENT:
            enabled = (True, False, False, False)
        elif t == self.MUSR_WAVELENGTH:
            enabled = (False, True, False, True)
        elif t == self.MUS_WAVELENGTH:
            enabled = (False, False, True, True)
        else:
            raise ValueError('Unexpected scattering value type!')

        self.solidContentParameter.setEnabled(enabled[0])
        self.musrParameter.setEnabled(enabled[1])
        self.musParameter.setEnabled(enabled[2])
        self.wavelengthParameter.setEnabled(enabled[3])

        self.definitionTypeComboBox.setCurrentIndex(t)

    def valueType(self) -> int:
        return self.definitionTypeComboBox.currentIndex()

    def setTolerancesVisible(self, visible: bool):
        self.solidContentParameter.setTolerancesVisible(visible)
        self.musrParameter.setTolerancesVisible(visible)
        self.musParameter.setTolerancesVisible(visible)

    def validateParticleRi(self, ri: RefractiveIndex, silent: bool = False) -> bool:
        if self.valueType() in (self.MUS_WAVELENGTH, self.MUSR_WAVELENGTH):
            if isinstance(ri, RefractiveIndex):
                if not ri.is_valid_wavelength(self.sciWavelength()):
                    if not silent:
                        low, high = ri.wrange
                        self.errorDialog.showMessage(
                            '{}:Wavelength: Refractive index of particles '
                            '(from {} to {} nm) is not '
                            'defined at the given wavelength!'.format(
                                self.groupBox.title(), low*1e9, high*1e9))
                        self.wavelengthParameter.focusValue()
                    return False
        return True

    def validateMediumRi(self, ri: RefractiveIndex, silent: bool = False) -> bool:
        if self.valueType() in (self.MUS_WAVELENGTH, self.MUSR_WAVELENGTH):
            if isinstance(ri, RefractiveIndex):
                if not ri.is_valid_wavelength(self.sciWavelength()):
                    if not silent:
                        low, high = ri.wrange
                        self.errorDialog.showMessage(
                            '{}:Wavelength: Refractive index of medium '
                            '(from {} to {} nm) is not '
                            'defined at the given wavelength!'.format(
                                self.groupBox.title(), low*1e9, high*1e9))
                        self.wavelengthParameter.focusValue()
                    return False
        return True

    def validate(self, prefix: str = None, silent: bool = False):
        if prefix is None:
            prefix = self.groupBox.title()
        vt = self.valueType()
        if vt == self.SOLID_CONTENT:
            if not self.validateSolidContent(prefix, silent=silent):
                return False
        elif vt == self.MUS_WAVELENGTH:
            if not self.validateMus(prefix, silent=silent):
                return False
            if not self.validateWavelength(prefix, silent=silent):
                return False
        elif vt == self.MUSR_WAVELENGTH:
            if not self.validateMusr(prefix, silent=silent):
                return False
            if not self.validateWavelength(prefix, silent=silent):
                return False

        return True

    def todict(self) -> dict:
        return {
            'type': self.__class__.__name__,
            'title': self.title(),
            'valueType': self.valueType(),
            'mus': self.musParameter.todict(),
            'musr': self.musrParameter.todict(),
            'solidContent': self.solidContentParameter.todict(),
            'wavelength': self.wavelengthParameter.todict(),
        }

    def fromdict(self, data: dict):
        t = data.get('type')
        if t != self.__class__.__name__:
            raise TypeError(
                'Expected data of type {} but got {}!'.format(
                    self.__class__.__name__, t))
        self.setTitle(data.get('title'))
        self.setValueType(data.get('valueType'))
        self.musParameter.fromdict(data.get('mus'))
        self.musrParameter.fromdict(data.get('musr'))
        self.solidContentParameter.fromdict(data.get('solidContent'))
        self.wavelengthParameter.fromdict(data.get('wavelength'))


class BaseSuspensionWidget(QScrollArea):
    def __init__(self, parent=None):
        super().__init__(parent)

        self.errorDialog = ErrorMessage(self)
        self.errorDialog.setWindowModality(Qt.WindowModality.WindowModal)

        self.warningDialog = WarningConfirmMessage(self)
        self.warningDialog.setWindowModality(Qt.WindowModality.WindowModal)

        self.particleSizeWidget = ParticleSizeWidget('Particle size')
        self.particleDensityWidget = DensityWidget(PARTICLE_MATERIALS, 'Particle density')
        self.particleRiWidget = RefractiveIndexWidget(PARTICLE_MATERIALS, 'Pearticle refractive index')
        self.mediumDensityWidget = DensityWidget(MEDIUM_MATERIALS, 'Medium density')
        self.mediumRiWidget = RefractiveIndexWidget(MEDIUM_MATERIALS, 'Medium refractive index')
        self.scatteringWidget = ScatteringWidget('Scattering properties')
        self.volumeWidget = NormalParameterWidget('Available volume', 'ml', 1e-6)
        self.volumeWidget.setTolerancesVisible(False)

        self.scatteringWidget.setValueType(self.scatteringWidget.SOLID_CONTENT)

        layout = QVBoxLayout()
        layout.addWidget(self.particleSizeWidget)
        layout.addWidget(self.particleDensityWidget)
        layout.addWidget(self.particleRiWidget)
        layout.addWidget(self.mediumDensityWidget)
        layout.addWidget(self.mediumRiWidget)
        layout.addWidget(self.scatteringWidget)
        layout.addWidget(self.volumeWidget)
        layout.addStretch()

        self.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)

        widget = QWidget()
        widget.setLayout(layout)
        self.setWidget(widget)
        self.setWidgetResizable(True)

        #self.setLayout(layout)

    def suspension(self) -> Suspension:
        pd = Normal(self.particleSizeWidget.sciDiameter(),
                    self.particleSizeWidget.sciDiameterStd())
        particle_ri = self.particleRiWidget.sciRefractiveIndex()
        medium_ri = self.mediumRiWidget.sciRefractiveIndex()
        particle_density = self.particleDensityWidget.sciDensity()
        susp = Suspension(pd, particle_ri=particle_ri, medium_ri=medium_ri,
                          particle_density=particle_density)
        t = self.scatteringWidget.valueType()
        if t == self.scatteringWidget.SOLID_CONTENT:
            susp.set_solid_content(self.scatteringWidget.sciSolidContent())
        elif t == self.scatteringWidget.MUSR_WAVELENGTH:
            susp.set_musr(self.scatteringWidget.musr(),
                          self.scatteringWidget.sciWavelength())
        elif t == self.scatteringWidget.MUS_WAVELENGTH:
            susp.set_mus(self.scatteringWidget.mus(),
                         self.scatteringWidget.sciWavelength())
        else:
            self.errorDialog.showMessage(
                'Unsupported definition of suspension scattering!')

        return susp

    def validate(self, silent: bool = False) -> bool:
        if not self.particleSizeWidget.validate(silent=silent):
            if not silent:
                self.ensureWidgetVisible(self.particleSizeWidget)
            return False
        if not self.particleDensityWidget.validate(silent=silent):
            if not silent:
                self.ensureWidgetVisible(self.particleDensityWidget)
            return False
        if not self.particleRiWidget.validate(silent=silent):
            if not silent:
                self.ensureWidgetVisible(self.particleRiWidget)
            return False
        if not self.mediumRiWidget.validate(silent=silent):
            if not silent:
                self.ensureWidgetVisible(self.mediumRiWidget)
            return False
        if not self.scatteringWidget.validate(silent=silent):
            if not silent:
                self.ensureWidgetVisible(self.scatteringWidget)
            return False
        if not self.volumeWidget.validate(silent=silent):
            if not silent:
                self.ensureWidgetVisible(self.volumeWidget)
            return False

        ri_particle = self.particleRiWidget.refractiveIndex()
        if not self.scatteringWidget.validateParticleRi(ri_particle, silent=silent):
            if not silent:
                self.ensureWidgetVisible(self.scatteringWidget)
            return False

        ri_medium = self.mediumRiWidget.refractiveIndex()
        if not self.scatteringWidget.validateParticleRi(ri_medium, silent=silent):
            if not silent:
                self.ensureWidgetVisible(self.scatteringWidget)
            return False

        m1 = self.particleRiWidget.material()
        m2 = self.particleDensityWidget.material()
        if m1 != MATERIAL_CUSTOM and m2 != MATERIAL_CUSTOM and m1 != m2:
            if not silent:
                result = self.warningDialog.showMessage(
                        'Different standard materials are used for particle '
                        'density ({}) and refractive index ({})! '
                        'Are the data corret?'.format(m2.name, m1.name))
                if result == QMessageBox.StandardButton.No:
                    return False
    
        m1 = self.mediumRiWidget.material()
        m2 = self.mediumDensityWidget.material()
        if m1 != MATERIAL_CUSTOM and m2 != MATERIAL_CUSTOM and m1 != m2:
            if not silent:
                result = self.warningDialog.showMessage(
                        'Different standard materials are used for medium '
                        'density ({}) and refractive index ({})! '
                        'Are the data corret?'.format(m2.name, m1.name))
                if result == QMessageBox.StandardButton.No:
                    return False

        return True

    def strDescriptor(self):
        d = self.particleSizeWidget.sciDiameter()
        d_std = self.particleSizeWidget.sciDiameterStd()
        particle = self.particleRiWidget.material()
        medium = self.mediumRiWidget.material()
        vt = self.scatteringWidget.valueType()
        if vt == self.scatteringWidget.SOLID_CONTENT:
            scattering = '{} {}'.format(self.scatteringWidget.sciSolidContent(), 'wt %')
        elif vt == self.scatteringWidget.MUSR_WAVELENGTH:
            musr = self.scatteringWidget.sciMusr()
            w = self.scatteringWidget.sciWavelength()
            scattering = "μs' {:.2f} 1/cm @ {:.0f} nm".format(musr*1e-2, w*1e9)
        elif vt == self.scatteringWidget.MUS_WAVELENGTH:
            mus = self.scatteringWidget.sciMus()
            w = self.scatteringWidget.sciWavelength()
            scattering = 'μs {:.2f} 1/cm @ {:.0f} nm'.format(mus*1e-2, w*1e9)

        volume = self.volumeWidget.sciValue()
        return '{:s} in {:s}, d {:.3f} ({:.3f}) μm, {:s}, V {:.1f} ml'.format(
            particle.name, medium.name,
            d*1e6, d_std*1e6, scattering, volume*1e6)

    def todict(self) -> dict:
        return {
            'type': self.__class__.__name__,

            'particleSize': self.particleSizeWidget.todict(),
            'particleDensity': self.particleDensityWidget.todict(),
            'particleRi': self.particleRiWidget.todict(),

            'mediumDensity': self.mediumDensityWidget.todict(),
            'mediumRi': self.mediumRiWidget.todict(),

            'scattering': self.scatteringWidget.todict(),
            'volume': self.volumeWidget.todict()
        }

    def fromdict(self, data: dict):
        t = data.get('type')
        if t != self.__class__.__name__:
            raise TypeError(
                'Expected data of type {} but got {}!'.format(
                    self.__class__.__name__, t))

        self.particleSizeWidget.fromdict(data.get('particleSize'))
        self.particleDensityWidget.fromdict(data.get('particleDensity'))
        self.particleRiWidget.fromdict(data.get('particleRi'))

        self.mediumDensityWidget.fromdict(data.get('mediumDensity'))
        self.mediumRiWidget.fromdict(data.get('mediumRi'))

        self.scatteringWidget.fromdict(data.get('scattering'))
        self.volumeWidget.fromdict(data.get('volume'))

    def sizeHint(self):
        return QSize(480, 640)

class DilutedSuspensionWidget(QScrollArea):
    def __init__(self, parent=None):
        super().__init__(parent)

        self.errorDialog = ErrorMessage(self)
        self.errorDialog.setModal(True)

        self.scatteringWidget = ScatteringWidget('Scattering properties')
        self.scatteringWidget.setTolerancesVisible(False)

        self.scatteringWidget.setValueType(self.scatteringWidget.MUSR_WAVELENGTH)

        self.volumeWidget = NormalParameterWidget('Required volume', 'ml', 1e-6)
        self.volumeWidget.setTolerancesVisible(False)
        self.volumeWidget.setValueValidator(NonzeroPositiveValidator())

        layout = QVBoxLayout()
        layout.addWidget(self.scatteringWidget)
        layout.addWidget(self.volumeWidget)
        layout.addStretch()

        self.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)

        widget = QWidget()
        widget.setLayout(layout)
        self.setWidget(widget)
        self.setWidgetResizable(True)

        #self.setLayout(layout)

    def validate(self, silent: bool = False):
        if not self.scatteringWidget.validate(silent=silent):
            return False
        if not self.volumeWidget.validate(silent=silent):
            return False

        return True

    def validateMediumRi(self, ri: RefractiveIndex, silent: bool = False) -> bool:
        return self.scatteringWidget.validateMediumRi(ri, silent=silent)

    def validateParticleRi(self, ri: RefractiveIndex, silent: bool = False) -> bool:
        return self.scatteringWidget.validateParticleRi(ri, silent=silent)

    def todict(self) -> dict:
        return {
            'type': self.__class__.__name__,
            'scattring': self.scatteringWidget.todict(),
            'volume': self.volumeWidget.todict()
        }

    def update(self, data: dict):
        self.scatteringWidget.update(data.get('scattering'))
        self.volumeWidget.update(data.get('volume'))

    def strDescriptor(self):
        vt = self.scatteringWidget.valueType()
        if vt == self.scatteringWidget.SOLID_CONTENT:
            scattering = '{} {}'.format(self.scatteringWidget.sciSolidContent(), 'wt %')
        elif vt == self.scatteringWidget.MUSR_WAVELENGTH:
            musr = self.scatteringWidget.sciMusr()
            w = self.scatteringWidget.sciWavelength()
            scattering = "μs' {:.2f} 1/cm @ {:.0f} nm".format(musr*1e-2, w*1e9)
        elif vt == self.scatteringWidget.MUS_WAVELENGTH:
            mus = self.scatteringWidget.sciMus()
            w = self.scatteringWidget.sciWavelength()
            scattering = 'μs {:.2f} 1/cm @ {:.0f} nm'.format(mus*1e-2, w*1e9)

        volume = self.volumeWidget.sciValue()
        return '{:s}, V {:.1f} ml'.format(scattering, volume*1e6)

    def todict(self) -> dict:
        return {
            'type': self.__class__.__name__,
            'scattering': self.scatteringWidget.todict(),
            'volume': self.volumeWidget.todict()
        }

    def updateSuspension(self, suspension: Suspension) -> Suspension:
        vt = self.scatteringWidget.valueType()
        if vt == self.scatteringWidget.SOLID_CONTENT:
            suspension.set_solid_content(self.scatteringWidget.sciSolidContent())
        elif vt == self.scatteringWidget.MUSR_WAVELENGTH:
            suspension.set_musr(self.scatteringWidget.sciMusr(),
                                self.scatteringWidget.sciWavelength())
        elif vt == self.scatteringWidget.MUSR_WAVELENGTH:
            suspension.set_musr(self.scatteringWidget.sciMus(),
                                self.scatteringWidget.sciWavelength())

        return suspension

    def fromdict(self, data: dict):
        t = data.get('type')
        if t != self.__class__.__name__:
            raise TypeError(
                'Expected data of type {} but got {}!'.format(
                    self.__class__.__name__, t))

        self.scatteringWidget.fromdict(data.get('scattering'))
        self.volumeWidget.fromdict(data.get('volume'))

    def sizeHint(self):
        return QSize(480, 640)


class ConfigurationDialog(QDialog):
    def __init__(self, widget: BaseSuspensionWidget or 
                               DilutedSuspensionWidget or None = None,
                 parent: QWidget = None):
        super().__init__(parent=parent)
        self.setModal(True)

        main_layout = QVBoxLayout()

        self.widgetLayout = QVBoxLayout()
        self.widgetLayout.setContentsMargins(0, 0, 0, 0)
        main_layout.addLayout(self.widgetLayout)

        self.okButton = QPushButton('Ok')
        self.okButton.clicked.connect(self.onOkButton)
        self.cancelButton = QPushButton('Cancel')
        self.cancelButton.clicked.connect(self.onCancelButton)

        buttonLayout = QHBoxLayout()
        buttonLayout.addStretch(1)
        buttonLayout.addWidget(self.okButton)
        buttonLayout.addWidget(self.cancelButton)

        main_layout.addLayout(buttonLayout)

        if widget:
            self.setWidget(widget)

        self.setLayout(main_layout)

    def setWidget(self, widget: BaseSuspensionWidget or
                                DilutedSuspensionWidget or None):
        if self.widgetLayout.count() > 0:
            item = self.widgetLayout.itemAt(0)
            currentWidget = item.widget()
            self.widgetLayout.removeWidget(currentWidget)
            currentWidget.setParent(None)

        if widget is not None:
            widget.setContentsMargins(0, 0, 0, 0)
            self.widgetLayout.addWidget(widget, 1)

    def widget(self) -> QWidget or None:
        if self.widgetLayout.count() > 0:
            return self.widgetLayout.itemAt(0).widget()

    def onOkButton(self):
        widget = self.widget()
        if widget is not None and widget.validate():
            self.accept()

    def onCancelButton(self):
        self.reject()


class DilutedSuspensionListItemWidget(QWidget):
    remove = QSignal(object)
    edit = QSignal(object)

    def __init__(self, listWidgetItem: QListWidgetItem,
                 suspension: DilutedSuspensionWidget):
        super().__init__()
        self.label = QLabel(suspension.strDescriptor())
        self._suspension = suspension
        self._listWidgetItem = listWidgetItem

        self.removeButton = QPushButton('Delete')
        self.removeButton.setIcon(
            self.removeButton.style().standardIcon(QStyle.SP_TrashIcon))
        self.removeButton.clicked.connect(lambda: self.remove.emit(self))

        self.editButton = QPushButton('Edit')
        self.editButton.setIcon(
            self.editButton.style().standardIcon(QStyle.SP_ArrowRight))
        self.editButton.clicked.connect(lambda: self.edit.emit(self))

        layout = QHBoxLayout()
        layout.addWidget(self.editButton)
        layout.addWidget(self.label)
        layout.addStretch(1)
        layout.addWidget(self.removeButton)

        self.setLayout(layout)

    def updateDescriptor(self):
        descriptor = self._suspension.strDescriptor()
        row = self.row()
        self.label.setText('({}) {}'.format(row + 1, descriptor))

    def listWidgetItem(self) -> QListWidgetItem:
        return self._listWidgetItem

    def setLisWidgetItem(self, lwi: QListWidgetItem or None):
        self._listWidgetItem = lwi

    def row(self) -> int or None:
        if self._listWidgetItem:
            lw = self._listWidgetItem.listWidget()
            return lw.row(self._listWidgetItem)

    def suspensionWidget(self) -> DilutedSuspensionWidget:
        return self._suspension


class DilutedSuspensionList(QListWidget):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.editSuspensionDialog = ConfigurationDialog(parent=self)

    def addSuspension(self, suspension: DilutedSuspensionWidget):
        listItem = QListWidgetItem()
        itemWidget = DilutedSuspensionListItemWidget(listItem, suspension)
        itemWidget.layout().setContentsMargins(1, 1, 1, 1)
        listItem.setSizeHint(itemWidget.sizeHint())
        suspension.listItem = listItem
        self.addItem(listItem)
        self.setItemWidget(listItem, itemWidget)
        itemWidget.remove.connect(
            lambda w: self.removeSuspension(w.suspensionWidget()))
        itemWidget.edit.connect(
            lambda w: self.editSuspension(w.suspensionWidget()))
        itemWidget.updateDescriptor()

    def removeSuspension(self, suspension: DilutedSuspensionWidget):
        listItem = getattr(suspension, 'listItem')
        if listItem:
            widget = self.itemWidget(listItem)
            widget.setLisWidgetItem(None)
            row = self.row(listItem)
            self.takeItem(row)
            suspension.listItem = None
            self.renumber(row)

    def editSuspension(self, suspension: DilutedSuspensionWidget):
        listItem = getattr(suspension, 'listItem')
        if listItem:
            self.editSuspensionDialog.setWidget(suspension)
            result = self.editSuspensionDialog.exec()
            self.editSuspensionDialog.setWidget(None)
            if not result and not suspension.validate(silent=True):
                self.takeItem(self.row(listItem))
                suspension.listItem = None
            else:
                self.itemWidget(listItem).updateDescriptor()

    def renumber(self, start: int = 0):
        for row in range(start, self.count()):
            widget = self.itemWidget(self.item(row))
            widget.updateDescriptor()
    
    def suspensions(self) -> Tuple[DilutedSuspensionWidget, ...]:
        return tuple(self.itemWidget(self.item(row)).suspensionWidget()
                     for row in range(self.count()))


def horizontalLine():
    hline = QFrame()
    hline.setFrameShape(QFrame.HLine)
    hline.setFrameShadow(QFrame.Sunken)
    return hline


class SequentialSuspensions(QWidget):
    def __init__(self):
        '''
        Widget that manages data of a sequential suspension dilution process.
        '''
        super().__init__()

        self.errorDialog = ErrorMessage(self)
        self.errorDialog.setModal(True)

        self.baseSuspensionWidget = BaseSuspensionWidget()
        self.baseSuspensionEditDialog = ConfigurationDialog(parent=self)

        self.baseSuspensionEditButton = QPushButton(
            'Create/Update base suspension')
        self.baseSuspensionEditButton.clicked.connect(
            self._updateBaseSuspension)
        self.baseSuspensionDescriptor = QLabel()

        self.addDilutedSuspensionButton = QPushButton(
            'Add a new diluted suspension')
        self.addDilutedSuspensionButton.clicked.connect(
            self._addDilutedSuspension)

        self.computeRecipeButton = QPushButton('Compute/update recipe')
        self.computeRecipeButton.clicked.connect(self.computeRecipe)

        self.dilutedSuspensionListWidget = DilutedSuspensionList()

        self.dilutionFeedbackMassRadioButton = QRadioButton(
            'Gravimetric feedback (g)')
        self.dilutionFeedbackVolumeRadioButton = QRadioButton(
            'Volumetric feedback (ml)')
        radio_layout = QHBoxLayout()
        radio_layout.setContentsMargins(0, 0, 0, 0)
        radio_layout.addStretch()
        radio_layout.addWidget(self.dilutionFeedbackMassRadioButton)
        radio_layout.addStretch()
        radio_layout.addWidget(self.dilutionFeedbackVolumeRadioButton)
        radio_layout.addStretch()
        self.dilutionFeedbackGroupBox = QGroupBox()
        self.dilutionFeedbackGroupBox.setTitle(
            'Select units of feedback data (blue && red columns)')
        self.dilutionFeedbackGroupBox.setLayout(radio_layout)
        self.dilutionFeedbackMassRadioButton.setChecked(True)
        self.dilutionFeedbackMassRadioButton.clicked.connect(
            lambda: self.setDilutionFededbackType('mass'))
        self.dilutionFeedbackVolumeRadioButton.clicked.connect(
            lambda: self.setDilutionFededbackType('volume'))

        self.dilutionTableWidget = QTableWidget(0, 8)
        self.dilutionTableWidget.setHorizontalHeaderLabels(
            ('Diluted Suspension',
             'Take (ml (g))', 'Taken (g)', 'Tol. (g)', 
             'Dilute (ml (g))', 'Diluted (g)', 'Tol. (g)',
             'Realization'))
        self.dilutionTableWidget.setSizeAdjustPolicy(
            QTableWidget.AdjustToContents)
        headerView = self.dilutionTableWidget.horizontalHeader()
        headerView.setSectionResizeMode(0, QHeaderView.ResizeToContents)
        headerView.setSectionResizeMode(1, QHeaderView.ResizeToContents)
        headerView.setSectionResizeMode(2, QHeaderView.ResizeToContents)
        headerView.setSectionResizeMode(3, QHeaderView.ResizeToContents)
        headerView.setSectionResizeMode(4, QHeaderView.ResizeToContents)
        headerView.setSectionResizeMode(5, QHeaderView.ResizeToContents)
        headerView.setSectionResizeMode(6, QHeaderView.ResizeToContents)
        headerView.setSectionResizeMode(7, QHeaderView.Stretch)
        headerView.setDefaultAlignment(Qt.AlignCenter)
        #.setSectionResizeMode(QHeaderView.ResizeToContents)

        self.analyzeRecipeFeedbackButton = QPushButton(
            'Analyze recipe feedback and compute realization')
        self.analyzeRecipeFeedbackButton.clicked.connect(
            self.analyzeRecipeFeedback)

        self.loadButton = QPushButton('Load')
        self.loadButton.setIcon(
            self.loadButton.style().standardIcon(QStyle.SP_DialogOpenButton))
        self.saveButton = QPushButton('Save')
        self.saveButton.setIcon(
            self.saveButton.style().standardIcon(QStyle.SP_DialogSaveButton))
        io_layout = QHBoxLayout()
        io_layout.setContentsMargins(0, 0, 0, 0)
        io_layout.addWidget(self.loadButton)
        io_layout.addStretch(1)
        io_layout.addWidget(self.saveButton)
        self.loadButton.clicked.connect(lambda x: self.loadFromFile())
        self.saveButton.clicked.connect(lambda x: self.saveToFile())

        layout = QVBoxLayout()

        gb = QGroupBox()
        gb_layout = QVBoxLayout()
        gb_layout.addWidget(QLabel('<h3><b>Base/given suspension</b></h3>'), 0, Qt.AlignCenter)
        gb_layout.addWidget(self.baseSuspensionEditButton)
        gb_layout.addWidget(self.baseSuspensionDescriptor, 0, Qt.AlignCenter)
        gb.setLayout(gb_layout)
        layout.addWidget(gb)

        #layout.addWidget(horizontalLine())
        #layout.addWidget(horizontalLine())
        gb = QGroupBox()
        gb_layout = QVBoxLayout()
        gb_layout.addWidget(QLabel('<h3><b>List of diluted suspensions</b></h3>'), 0, Qt.AlignCenter)
        gb_layout.addWidget(self.addDilutedSuspensionButton)
        gb_layout.addWidget(self.dilutedSuspensionListWidget)
        gb.setLayout(gb_layout)
        layout.addWidget(gb)

        #layout.addWidget(horizontalLine())
        #layout.addWidget(horizontalLine())
        gb = QGroupBox()
        gb_layout = QVBoxLayout()
        gb_layout.addWidget(QLabel('<h3><b>Sequential dilution recipe</b></h3>'), 0, Qt.AlignCenter)
        gb_layout.addWidget(self.computeRecipeButton)
        # layout.addWidget(QLabel('<b>Select units of feedback data (blue && red columns)</b>'), 0, Qt.AlignCenter)
        gb_layout.addWidget(self.dilutionFeedbackGroupBox)
        gb_layout.addWidget(self.dilutionTableWidget)
        gb_layout.addWidget(self.analyzeRecipeFeedbackButton)
        gb.setLayout(gb_layout)
        layout.addWidget(gb)

        # layout.addWidget(horizontalLine())
        # layout.addWidget(horizontalLine())
        layout.addLayout(io_layout)
    
        self.setLayout(layout)

    def _updateBaseSuspension(self):
        '''
        Show a dialog for entering/updating the data of base suspension and
        process the dialog result.
        '''
        label = ''
        self.baseSuspensionEditDialog.setWidget(self.baseSuspensionWidget)
        result = self.baseSuspensionEditDialog.exec()
        self.baseSuspensionEditDialog.setWidget(None)
        if result or self.baseSuspensionWidget.validate(silent=True):
            label = self.baseSuspensionWidget.strDescriptor()
        self.baseSuspensionDescriptor.setText(label)

    def _addDilutedSuspension(self):
        '''
        Show a dialog for entering a new diluted suspension and process the
        dialog result. 
        '''
        suspension = DilutedSuspensionWidget()
        dialog = ConfigurationDialog(suspension, self)
        dialog.setWindowTitle('Diluted suspension')
        if dialog.exec():
            suspension.setParent(None)
            self.dilutedSuspensionListWidget.addSuspension(suspension)

    def validateBaseSuspension(self, silent: bool = False) -> bool:
        '''
        Validate the data of the base suspension.

        Parameters
        ----------
        silent: bool
            Do not show error dialogs if True.

        Returns
        -------
        valid: bool
            True if the data of the base suspension are valid, else False.
        '''
        # validate base suspension
        if not self.baseSuspensionDescriptor.text():
            if not silent:
                self.errorDialog.showMessage('Base/given suspension is not defined!')
            return False

        if not self.baseSuspensionWidget.validate(silent=True):
            if not silent:
                self.errorDialog.showMessage('Base/given suspension is not valid!')
            return False

        return True

    def validateDilutedSuspensions(self, silent: bool = False) -> bool:
        '''
        Validate the data of all the listed diluted Suspensions.

        Parameters
        ----------
        silent: bool
            Do not show error dialogs if True.

        Returns
        -------
        valid: bool
            True if the date of diluted suspensions are valid, else False.

        Note
        ----
        Requires valid base suspension.
        '''
        if not self.validateBaseSuspension():
            return False

        baseSuspensionWidget = self.baseSuspensionWidget

        dilutedSuspensionWidgets = self.dilutedSuspensionListWidget.suspensions()
        if not dilutedSuspensionWidgets:
            if not silent:
                self.errorDialog.showMessage(
                    'Define at least one diluted suspension!')
            return False

        # validate the individual diluted suspensions against the base
        ri_particle = baseSuspensionWidget.particleRiWidget.refractiveIndex()
        ri_medium = baseSuspensionWidget.mediumRiWidget.refractiveIndex()
        for index, reqSuspWidget in enumerate(dilutedSuspensionWidgets):
            wavelength = reqSuspWidget.scatteringWidget.sciWavelength()
            valid = reqSuspWidget.validateParticleRi(ri_particle, silent=True)
            if not valid:
                if not silent:
                    w1, w2 = ri_particle.wrange
                    self.errorDialog.showMessage(
                        'Wavelength ({:.0f} nm) property of the diluted '
                        'suspension ({}) is not within the valid range '
                        '({:.0f} - {:.0f} nm) of the particle refractive '
                        'index!'.format(
                            wavelength*1e9, index + 1, w1*1e9, w2*1e9))
                return False
            valid = reqSuspWidget.validateMediumRi(ri_medium, silent=True)
            if not valid:
                if not silent:
                    w1, w2 = ri_particle.wrange
                    self.errorDialog.showMessage(
                        'Wavelength ({:.0f} nm) property of the diluted '
                        'suspension ({}) is not within the valid range '
                        '({:.0f} - {:.0f} nm) of the medium refractive '
                        'index!'.format(
                            wavelength*1e9, index + 1, w1*1e9, w2*1e9))
                return False

        return True

    def baseSuspension(self, silent: bool = False) \
            -> Tuple[Suspension, float] or None:
        '''
        Returns the base suspension as a :py:class:`Suspension` instances
        and the available volume of the base suspension (m3).

        Parameters
        -----------
        silent: bool
            Show no error dialogs if True.

        Returns
        -------
        suspension: Suspension
            Base suspension as a :py:class:`Suspension` instance.
        volume: float
            Available volume of the base suspension (m3).

        Note
        ----
        Returns a single None value if data of the base suspension are invalid
        or not available.
        '''
        result = None
        if self.validateBaseSuspension(silent=silent):
            suspension = self.baseSuspensionWidget.suspension()
            volume = self.baseSuspensionWidget.volumeWidget.sciValue()
            result = (suspension, volume)

        return result

    def dilutedSuspensions(self, silent: bool = False) \
            -> Tuple[Tuple[Suspension, ...], Tuple[float, ...]] or None:
        '''
        Generate and return a tuple of :py:class:`Suspension` instances and
        target volumnes for the current list of diluted suspensions.

        Returns
        -------
        suspensions: Tuple[Suspension, ...]
            A tuple of diluted suspensions that follows the order of
            suspensions in the suspension list (not dilution table).
        volumes: Tuple[float, ...]
            Corresponding target volumes of diluted suspensions (m3). 
        '''
        if not self.validateBaseSuspension(silent=silent):
            return None
        if not self.validateDilutedSuspensions(silent=silent):
            return None

        dilutedSuspensions = []
        dilutedVolumes = []
        dilutedSuspensionWidgets = self.dilutedSuspensionListWidget.suspensions()
        baseSuspension, _ = self.baseSuspension(silent=True)
        for index, dilutedSuspensionWidget in enumerate(dilutedSuspensionWidgets):
            suspension, _ = self.baseSuspension(silent=True)
            dilutedSuspensionWidget.updateSuspension(suspension)
            if suspension.solid_content() > baseSuspension.solid_content():
                if not silent:
                    self.errorDialog.showMessage(
                        'Solid content of the diluted suspension ({}) '
                        '({:.2f} wt %) exceeds the solid content '
                        '({:.2f} wt %) of the given suspension!'.format(
                            index + 1, suspension.solid_content(),
                            baseSuspension.solid_content()))
                    return None
            dilutedSuspensions.append(suspension)
            dilutedVolumes.append(dilutedSuspensionWidget.volumeWidget.sciValue())

        return (dilutedSuspensions, dilutedVolumes)

    def computeRecipe(self, silent: bool = False, feedback: bool = False):
        '''
        Compute dilution recipe and fill the dilution table. Opionally
        process feedback data and fill the data into the dilution table.

        Parameters
        ----------
        silent: bool
            If True, show no error dialogs, else show error dialogs.
        feedback: bool
            If True, process dilution feedback data.
        '''
        # validate the base suspension
        if not self.validateBaseSuspension():
            return None

        # validate the required/diluted suspensions
        if not self.validateDilutedSuspensions():
            return None

        base = self.baseSuspension(silent=silent)
        if base is None:
            return None
        baseSuspension, baseSuspensionVolume = base

        diluted = self.dilutedSuspensions(silent=silent)
        if diluted is None:
            return None
        dilutedSuspensions, dilutedSuspensionsVolume = diluted

        # compute dilution recipe
        volumeRecipe, massRecipe, indices = \
            baseSuspension.sequential_dilution_recipe(
                dilutedSuspensions, dilutedSuspensionsVolume)
        # sort diluted suspensions to match the recipe
        sortedSuspensions = tuple(dilutedSuspensions[indx] for indx in indices)

        requiredBaseSuspensionVolume = volumeRecipe[0][0]
        if baseSuspensionVolume < requiredBaseSuspensionVolume:
            if not silent:
                self.errorDialog.showMessage(
                    'Volume of the base suspenson ({:.3f} ml) is insufficient to '
                    'prepare the given dilutions ({:.3f} ml)!'.format(
                        baseSuspensionVolume*1e6,
                        requiredBaseSuspensionVolume*1e6))

        dilutedSuspensionWidgets = self.dilutedSuspensionListWidget.suspensions()
        tableWidget = self.dilutionTableWidget
        # populate the mix table with results
        self._prepareDilutionTable(tableWidget, len(dilutedSuspensions))
        for row, (vr, mr, s) in enumerate(zip(volumeRecipe, massRecipe, sortedSuspensions)):
            v_take, v_diluted = vr
            m_take, m_diluted = mr
            index = dilutedSuspensions.index(s)
            widget = dilutedSuspensionWidgets[index]
            # Suspension descriptor string
            item = tableWidget.item(row, 0)
            item.setText('({}) {}'.format(index + 1, widget.strDescriptor()))
            item.setData(Qt.UserRole, (widget, s))
            # Volume and mass to take from the previous suspension
            item = tableWidget.item(row, 1)
            item.setText('{:.4f} ({:.4f})'.format(v_take*1e6, m_take*1e3))
            item.setData(Qt.UserRole, (v_take, m_take))
            # Dilute the taken volume/mass to the specified volume/mass
            item = tableWidget.item(row, 4)
            item.setText('{:.4f} + {:.4f} = {:.4f} ({:.4f})'.format(
                v_take*1e6, (v_diluted - v_take)*1e6, v_diluted*1e6, m_diluted*1e3))
            item.setData(Qt.UserRole, (v_diluted, m_diluted))

        if feedback:
            feedbackData, errors = self.collectDilutionFeedback(silent=silent)
            if not errors:
                feedbackType = self.dilutionFeedbackType()
                recipe = tuple((taken[0], diluted[0])
                               for taken, diluted in feedbackData)
                if feedbackType == 'mass':
                    realizationSuspensions, realizationMasses = \
                        baseSuspension.apply_sequential_mass_dilution_recipe(recipe)
                else:
                    realizationSuspensions, realizationVolumes = \
                        baseSuspension.apply_sequential_volume_dilution_recipe(recipe)
                for row, realization in enumerate(realizationSuspensions):
                    widget, target = tableWidget.item(row, 0).data(Qt.UserRole)
                    data, text = None, ''
                    vt = widget.scatteringWidget.valueType()
                    if vt == widget.scatteringWidget.SOLID_CONTENT:
                        sc_target = target.solid_content()
                        sc = realization.solid_content()
                        err = (sc - sc_target)/sc_target
                        text = '{:.2f} wt % ({:.2f}%)'.format(sc, err*100.0)
                        data = (realization, 'sc', sc_target, sc, err)
                    elif vt == widget.scatteringWidget.MUS_WAVELENGTH:
                        wavelength = widget.scatteringWidget.sciWavelength()
                        mus_target = target.mus(wavelength)
                        mus = realization.mus(wavelength)
                        err = (mus - mus_target)/mus_target
                        text = 'μs {:.2f} 1/cm ({:.2f}%)'.format(mus*1e-2, err*100.0)
                        data = (realization, 'mus', mus_target, mus, err)
                    elif vt == widget.scatteringWidget.MUSR_WAVELENGTH:
                        wavelength = widget.scatteringWidget.sciWavelength()
                        musr_target = target.musr(wavelength)
                        musr = realization.musr(wavelength)
                        err = (musr - musr_target)/musr_target
                        text = "μs' {:.2f} 1/cm ({:.2f}%)".format(musr*1e-2, err*100.0)
                        data = (realization, 'musr', musr_target, musr, err)
                    item = tableWidget.item(row, 7)
                    item.setText(text)
                    item.setData(Qt.UserRole, data)

    def analyzeRecipeFeedback(self, silent: bool = False):
        self.computeRecipe(silent=silent, feedback=True)

    def _checkRecipeDilutionFeedbackComplete(
            self, feedback: Tuple[Tuple[Tuple[float, float],
                                        Tuple[float, float]], ...]) -> bool:
        for take, dilute in feedback:
            if take[0] is None or dilute[0] is None:
                return False
            if take[0] < 0 or dilute[0] < 0:
                return False
            if take[0] > dilute[0]:
                return False

        return True

    def _chckRecipeRecipeToleranceFeedbackComplete(
            self, feedback: Tuple[Tuple[Tuple[float, float],
                                  Tuple[float, float]], ...]) -> bool:
        for take, dilute in feedback:
            if take[1] is None or dilute[1] is None:
                return False
            if take[1] < 0 or dilute[1] < 0:
                return False

        return True

    def _prepareDilutionTable(self, table: QTableWidget, n: int):
        '''
        Initialize the dilution table
        '''
        table.setRowCount(n)
        for row in range(n):
            for col in range(table.columnCount()):
                item = table.item(row, col)
                if item is None:
                    item = QTableWidgetItem('')
                    if col in (3, 6):
                        item.setText('0.0000') # tolerances
                    table.setItem(row, col, item)
                    if col == 0:
                        item.setTextAlignment(Qt.AlignCenter | Qt.AlignVCenter)
                        item.setFlags(item.flags() ^ Qt.ItemFlag.ItemIsEditable)
                    elif col in (1, 4):
                        item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                        item.setBackground(QColor(0, 255, 0, 32))
                        item.setFlags(item.flags() ^ Qt.ItemFlag.ItemIsEditable)
                    elif col in (2, 5):
                        item.setTextAlignment(Qt.AlignLeft | Qt.AlignVCenter)
                        item.setBackground(QColor(0, 0, 255, 32))
                    elif col in (3, 6):
                        item.setTextAlignment(Qt.AlignLeft | Qt.AlignVCenter)
                        item.setBackground(QColor(255, 0, 0, 32))
                    elif col == 7:
                        item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                        item.setBackground(QColor(255, 128, 0, 32))
                        item.setFlags(item.flags() ^ Qt.ItemFlag.ItemIsEditable)

    def collectDilutionFeedback(self, silent: bool = False) -> \
            Tuple[
                Tuple[Tuple[Tuple[float or None, float or None],
                            Tuple[float or None, float or None]], ...],
                      bool]:
        '''
        Collect dilution feedback data from the dilution table.

        Returns
        -------
        data: Tuple[Tuple[Tuple[float or None, float or None],
                          Tuple[float or None, float or None]], ...]
            Dilution feedback data for each dilution as a tupe
            (((take, take_tol),(dilute, dilute_tol)), ...). The values are
            returned in ntive units, i.e. in mg if feedback type is
            gravimetric ("mass") or ml if feedback type is volumetric
            ("volume").
        error: bool
            True if feedback data are incomplete or contain errors.
        '''
        tableWidget = self.dilutionTableWidget
        rows = tableWidget.rowCount()
        feedback = []
        errors = []
        nonzeroPositiveValidator = NonzeroPositiveValidator()
        nonnegativeValidator = NonnegativeValidator()
        for row in range(rows):
            susp = tableWidget.item(row, 0).text()
            taken = taken_tol = diluted = diluted_tol = None
            try:
                taken = float(tableWidget.item(row, 2).text())
                nonzeroPositiveValidator(taken)
            except:
                errors.append('Taken feedback value of suspension '
                              '"{}" is not valid!'.format(susp))
            try:
                taken_tol = float(tableWidget.item(row, 3).text())
                nonnegativeValidator(taken_tol)
            except:
                errors.append('Taken feedback tolerance of suspension '
                              '"{}" is not valid!'.format(susp))
            try:
                diluted = float(tableWidget.item(row, 5).text())
                nonzeroPositiveValidator(diluted)
            except:
                errors.append('Diluted to feedback value of suspension '
                              '"{}" is not valid!'.format(susp))
            try:
                diluted_tol = float(tableWidget.item(row, 6).text())
                nonnegativeValidator(diluted_tol)
            except:
                errors.append('Diluted to feedback tolerance of suspension '
                              '"{}" is not valid!'.format(susp))

            if taken is not None and diluted is not None and taken > diluted:
                errors.append('Suspension "{}": The value of taken feedback '
                              'exceeds the value of diluted '
                              'to feedback!'.format(susp))

            feedback.append(
                ((taken, taken_tol), (diluted, diluted_tol))
            )
        if not silent and errors:
            self.errorDialog.showMessage('\n'.join(errors))

        return tuple(feedback), bool(errors)

    def setDilutionFededbackType(self, t: Literal['mass', 'volume']):
        '''
        Set mixing table feedback information type (mass or volume)
        '''
        
        if t == 'mass':
            self.dilutionFeedbackMassRadioButton.setChecked(True)
            taken, taken_tol = 'Taken (g)', 'Tol. (g)'
            diluted, diluted_tol = 'Diluted (g)', 'Tol. (g)'
        elif t == 'volume':
            self.dilutionFeedbackVolumeRadioButton.setChecked(True)
            taken, taken_tol = 'Taken (ml)', 'Tol. (ml)'
            diluted, diluted_tol = 'Diluted (ml)', 'Tol. (ml)'
        else:
            self.errorDialog.showMessage(
                'Unexpected dilution feedback type "{}"!'.format(t))
            return 

        self.dilutionTableWidget.horizontalHeaderItem(2).setText(taken)
        self.dilutionTableWidget.horizontalHeaderItem(3).setText(taken_tol)
        self.dilutionTableWidget.horizontalHeaderItem(5).setText(diluted)
        self.dilutionTableWidget.horizontalHeaderItem(6).setText(diluted_tol)

    def dilutionFeedbackType(self) -> Literal['mass', 'volume']:
        '''
        Returns
        -------
        feedbackType: Literal['mass', 'volume']
            Dilution feedback data type/units.
        '''
        if self.dilutionFeedbackMassRadioButton.isChecked():
            return 'mass'
        else:
            return 'volume'

    def todict(self) -> dict:
        '''
        Export widget data and state to a dict.

        Returns
        -------
        data : dict
            Exportet widget data and state.
        '''
        return {
            'dilutionFeedbackType': self.dilutionFeedbackType(),
            'baseSuspension': self.baseSuspensionWidget.todict(),
            'dilutedSuspensions': [susp.todict() for susp in
                                   self.dilutedSuspensionListWidget.suspensions()],
            'dilutionFeedback': self.collectDilutionFeedback(silent=True)[0]
        }

    def fromdict(self, data: dict):
        '''
        Load widget data and state from a dict that was produced
        :py:meth:`~SequentialSuspensions.todict`. 

        Parameters
        ----------
        data: dict
            Widget data and state to be loaded.
        '''
        self.dilutionTableWidget.clearContents()
        baseSuspension = data.get('baseSuspension')
        if baseSuspension is not None:
            self.baseSuspensionWidget.fromdict(baseSuspension)
            label = ''
            if self.baseSuspensionWidget.validate():
                label = self.baseSuspensionWidget.strDescriptor()
            self.baseSuspensionDescriptor.setText(label)
        self.dilutedSuspensionListWidget.clear()
        dilutedSuspensions = data.get('dilutedSuspensions')
        for dilutedSuspension in dilutedSuspensions:
            widget = DilutedSuspensionWidget()
            widget.fromdict(dilutedSuspension)
            self.dilutedSuspensionListWidget.addSuspension(widget)
        if dilutedSuspensions:
            self.computeRecipe()
        self.setDilutionFededbackType(data.get('dilutionFeedbackType'))
        dilutionFeedback = data.get('dilutionFeedback', 'mass')
        if dilutionFeedback:
            tableWidget = self.dilutionTableWidget
            rows = self.dilutionTableWidget.rowCount()
            for row, feedback in zip(range(rows), dilutionFeedback):
                (taken, taken_tol), (diluted, diluted_tol) = feedback
                if taken is not None:
                    tableWidget.item(row, 2).setText('{:.4f}'.format(taken))
                if taken_tol is not None:
                    tableWidget.item(row, 3).setText('{:.4f}'.format(taken_tol))
                if diluted is not None:
                    tableWidget.item(row, 5).setText('{:.4f}'.format(diluted))
                if diluted_tol is not None:
                    tableWidget.item(row, 6).setText('{:.4f}'.format(diluted_tol))

            if self._checkRecipeDilutionFeedbackComplete(dilutionFeedback):
                self.analyzeRecipeFeedback()

    def saveToFile(self, filename: str or None = None):
        '''
        Save widget data and state to a JSON file.

        Parameters
        ----------
        filename: str
            Target JSON file. A file save dialog is shown if filename is None.
        '''
        if filename is None:
            filename = QFileDialog.getSaveFileName(
                self, 'Save configuration to JSON file', None,
                'JSON (*.json)')
            if filename:
                filename = filename[0]

        if filename:
            try:
                with open(filename, 'w') as fid:
                    json.dump(self.todict(), fid, indent=4)
            except Exception:
                self.errorDialog.showMessage(
                    'Failed to save configuration to "{}"!'.format(filename),
                    traceback.format_exc())

    def loadFromFile(self, filename: str or None = None):
        '''
        Load widget data and state from a JSON file.

        Parameters
        ----------
        filename: str
            JSON file to load. A file open dialog is shown if filename is None.
        '''
        if filename is None:
            filename = QFileDialog.getOpenFileName(
                self, 'Save configuration to JSON file', None,
                'JSON (*.json)')
            if filename:
                filename = filename[0]

        if filename:
            try:
                with open(filename, 'r') as fid:
                    data = json.load(fid)
                    self.fromdict(data)
            except Exception:
                self.errorDialog.showMessage(
                    'Failed to load configuration from "{}"!'.format(filename),
                    traceback.format_exc())

    def sizeHint(self) -> QSize:
        return QSize(1024, 720)


class SequentialVolumeDilutionRecipe():
    def __init__(self, sequence: Tuple[Tuple[float or Callable[[], float],
                                             float or Callable[[], float]], ...]):
        '''
        Sequence of suspension dilutions defined as (take, dilute), where take
        is the volume (m3) taken from a suspension and dilute the diluted
        volune (m3) of the diluted suspension. The diluted suspension is
        used as the base for the next dilution.

        Parameters
        ----------
        sequence: Tuple[Tuple[float or Callable[[], float], float or Callable[[], float]], ...]
            Dilution sequence as a tuple of tuples (take, dilute), where take
            is the volume (m3) taken from a suspension and dilute
            the diluted volune (m3). If take or dilute are callabels, they
            should take no parameters and return the corresponding volumes.
        
        '''
        self._dilution_sequence = sequence

    def apply(self, suspension: Suspension):
        '''
        Apply the dilution recipe to the given suspension.

        Parameters
        ----------
        suspension: Suspension
            Initial suspension of the dilution sequence.

        Returns
        -------
        suspensions: Tuple[Suspension, ...]
            A tuple of diluted suspensions.
        '''
        suspensions = []
        current = suspension
        for v_take, v_diluted in self._dilution_sequence:
            if callable(v_take):
                v_take = v_take()
            if callable(v_diluted):
                v_diluted = v_diluted()
            current = current.dilute_volume(v_take, v_diluted)
            suspensions.append(current)

        return tuple(suspensions)


class SequentialMassDilutionRecipe():
    def __init__(self, sequence: Tuple[Tuple[float or Callable[[], float],
                                             float or Callable[[], float]], ...]):
        '''
        Sequence of suspension dilutions defined as (take, dilute), where take
        is the mass (kg) taken from a suspension and dilute the diluted
        mass (kg).

        Parameters
        ----------
        sequence: Tuple[Tuple[float, float], ...]
            Dilution sequence as a tuple of tuples (take, dilute), where take
            is the mass (kg) taken from a suspension and dilute
            the diluted mass (kg). If take or dilute are callabels, they
            should take no parameters and return the corresponding mass.
        
        '''
        self._dilution_sequence = sequence

    def apply(self, suspension: Suspension):
        '''
        Apply the dilution recipe to the given suspension.

        Parameters
        ----------
        suspension: Suspension
            Initial suspension of the dilution sequence.

        Returns
        -------
        suspensions: Tuple[Suspension, ...]
            A tuple of diluted suspensions.
        '''
        suspensions = []
        current = suspension
        for m_take, m_diluted in self._dilution_sequence:
            if callable(m_take):
                m_take = m_take()
            if callable(m_diluted):
                m_diluted = m_diluted()
            current = current.dilute_volume(m_take, m_diluted)
            suspensions.append(current)

        return tuple(suspensions)

if __name__ == '__main__':
    import sys
    import os.path
    from xopto import DATA_PATH

    app = QApplication(sys.argv)

    icon = QIcon(os.path.join(DATA_PATH, 'xopto', 'xopto_square_bright.png'))
    app.setWindowIcon(icon)

    # ri = RefractiveIndexWidget()
    # ri.show()

    #ps = ParticleSizeWidget()
    #ps.show()

    #s = ScatteringWidget()
    #s.show()

    #sc = QScrollArea()
    #w = BaseSuspensionWidget()
    #sc.setWidget(gs)
    #sc.show()

    #w = DilutedSuspensionWidget()
    #w = QScrollArea()
    #w = BaseSuspensionWidget()
    #w.setWidget(rs)

    #w = ConfigurationDialog(w)
    #w.show()

    w = SequentialSuspensions()
    w.setWindowTitle('Sequential Suspension Dilution Data (SSDD)')
    w.show()

    #lw = QListWidget()
    #item = QListWidgetItem()
    #lw.addItem(item)
    #lw.setItemWidget(item, QLabel('Shit'))
    #lw.show()

    app.exec()

