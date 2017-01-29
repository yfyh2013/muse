"""
  MUSE -- A Multi-algorithm-collaborative Universal Structure-prediction Environment

  Copyright (C) 2010-2017 by Zhong-Li Liu

  This program is free software; you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software Foundation
  version 2 of the License.

  This program is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
  PARTICULAR PURPOSE.  See the GNU General Public License for more details.

  E-mail: zl.liu@163.com
"""

import math,warnings
import numpy as np
from numpy import dot
from numpy.linalg import norm

AtomicNumbers = {}
ChemicalSymbols = ['X',  'H',  'He', 'Li', 'Be', 'B',  'C',  'N',  'O',  'F',
                    'Ne', 'Na', 'Mg', 'Al', 'Si', 'P',  'S',  'Cl', 'Ar', 'K',
                    'Ca', 'Sc', 'Ti', 'V',  'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu',
                    'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y',
                    'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In',
                    'Sn', 'Sb', 'Te', 'I',  'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr',
                    'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm',
                    'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au',
                    'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac',
                    'Th', 'Pa', 'U',  'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es',
                    'Fm', 'Md', 'No', 'Lr']

for a, symbol in enumerate(ChemicalSymbols):
    AtomicNumbers[symbol] = a

names = {'position': ('positions', np.zeros(3)),'number':   ('numbers',   0)}

def atomproperty(name, doc):
    def getter(self):
        return self.get(name)

    def setter(self, value):
        self.set(name, value)

    def deleter(self):
        self.delete(name)

    return property(getter, setter, deleter, doc)


def xyzproperty(index):
    def getter(self):
        return self.position[index]

    def setter(self, value):
        self.position[index] = value

    return property(getter, setter, doc='XYZ'[index] + '-coordinate')


class Treat_Atom(object):
    __slots__ = ['data', 'atoms', 'index']

    def __init__(self, symbol='X', position=(0, 0, 0),atoms=None,index=None):

        self.data = d = {}

        if atoms is None:
            if isinstance(symbol, str):
                d['number'] = AtomicNumbers[symbol]
            else:
                d['number'] = symbol
            d['position'] = np.array(position, float)

        self.index = index
        self.atoms = atoms

    def __repr__(self):
        s = "Atom('%s', %s" % (self.symbol, list(self.position))

        if self.atoms is None:
            s += ')'
        else:
            s += ', index=%d)' % self.index
        return s

    def cut_reference_to_atoms(self):
        for name in names:
            self.data[name] = self.get_raw(name)
        self.index = None
        self.atoms = None

    def get_raw(self, name):
        if name == 'symbol':
            return ChemicalSymbols[self.get_raw('number')]

        if self.atoms is None:
            return self.data[name]

        plural = names[name][0]
        if plural in self.atoms.arrays:
            return self.atoms.arrays[plural][self.index]
        else:
            return None

    def get(self, name):
        value = self.get_raw(name)
        if value is None:
            if name == 'mass':
                value = atomic_masses[self.number]
            else:
                value = names[name][1]
        return value

    def set(self, name, value):
        if name == 'symbol':
            name = 'number'
            value = AtomicNumbers[value]

    def delete(self, name):
        assert self.atoms is None
        assert name not in ['number', 'symbol', 'position']
        self.data[name] = None

    symbol = atomproperty('symbol', 'Chemical symbol')
    number = atomproperty('number', 'Atomic number')
    position = atomproperty('position', 'XYZ-coordinates')

    x = xyzproperty(0)
    y = xyzproperty(1)
    z = xyzproperty(2)

    def _get(self, name):
        warnings.warn('Use atom.%s' % name, stacklevel=3)
        return getattr(self, name)

    def _set(self, name, value):
        warnings.warn('Use atom.%s = ...' % name, stacklevel=3)
        setattr(self, name, value)

    def get_symbol(self): return self._get('symbol')
    def get_atomic_number(self): return self._get('number')
    def get_position(self): return self._get('position')

    def set_symbol(self, value): self._set('symbol', value)
    def set_atomic_number(self, value): self._set('number', value)
    def set_position(self, value): self._set('position', value)



def unit_vector(x):
    y = np.array(x, dtype='float')
    return y/norm(y)

def cell_to_cellpar(cell):
    va, vb, vc = cell
    a = np.linalg.norm(va)
    b = np.linalg.norm(vb)
    c = np.linalg.norm(vc)
    alpha = 180.0/math.pi*np.arccos(dot(vb, vc)/(b*c))
    beta  = 180.0/math.pi*np.arccos(dot(vc, va)/(c*a))
    gamma = 180.0/math.pi*np.arccos(dot(va, vb)/(a*b))
    return np.array([a, b, c, alpha, beta, gamma])


def cellpar_to_cell(cellpar, ab_normal=(0,0,1), a_direction=None):
    if a_direction is None:
        if np.linalg.norm(np.cross(ab_normal, (1,0,0))) < 1e-5:
            a_direction = (0,0,1)
        else:
            a_direction = (1,0,0)

    ad = np.array(a_direction)
    Z = unit_vector(ab_normal)
    X = unit_vector(ad - dot(ad, Z)*Z)
    Y = np.cross(Z, X)

    alpha, beta, gamma = 90., 90., 90.
    if isinstance(cellpar, (int, long, float)):
        a = b = c = cellpar
    elif len(cellpar) == 1:
        a = b = c = cellpar[0]
    elif len(cellpar) == 3:
        a, b, c = cellpar
        alpha, beta, gamma = 90., 90., 90.
    else:
        a, b, c, alpha, beta, gamma = cellpar
    alpha *= math.pi/180.0
    beta *= math.pi/180.0
    gamma *= math.pi/180.0
    va = a * np.array([1, 0, 0])
    vb = b * np.array([math.cos(gamma), math.sin(gamma), 0])
    cx = math.cos(beta)
    cy = (math.cos(alpha) - math.cos(beta)*math.cos(gamma))/math.sin(gamma)
    #print cellpar
    #print 1. - cx*cx - cy*cy
    cz = math.sqrt(1. - cx*cx - cy*cy)

    vc = c * np.array([cx, cy, cz])

    abc = np.vstack((va, vb, vc))
    T = np.vstack((X, Y, Z))
    cell = dot(abc, T)
    return cell

class build_cell(object):
    def __init__(self, symbols=None,
                 positions=None, numbers=None,
                 tags=None, masses=None,
                 scaled_positions=None,
                 cell=None, pbc=None):

        atoms = None

        if hasattr(symbols, 'get_positions'):
            atoms = symbols
            symbols = None
        elif (isinstance(symbols, (list, tuple)) and
              len(symbols) > 0 and isinstance(symbols[0], Treat_Atom)):
            data = [[atom.get_raw(name) for atom in symbols]
                    for name in
                    ['position', 'number', 'tag', 'momentum',
                     'mass', 'magmom', 'charge']]
            atoms = self.__class__(None, *data)
            symbols = None

        if atoms is not None:
            if scaled_positions is not None:
                raise NotImplementedError
            if symbols is None and numbers is None:
                numbers = atoms.get_atomic_numbers()
            if positions is None:
                positions = atoms.get_positions()
            if tags is None and atoms.has('tags'):
                tags = atoms.get_tags()
            if cell is None:
                cell = atoms.get_cell()
            if pbc is None:
                pbc = atoms.get_pbc()

        self.arrays = {}

        if symbols is None:
            if numbers is None:
                if positions is not None:
                    natoms = len(positions)
                elif scaled_positions is not None:
                    natoms = len(scaled_positions)
                else:
                    natoms = 0
                numbers = np.zeros(natoms, int)
            self.new_array('numbers', numbers, int)
        else:
            if numbers is not None:
                raise ValueError(
                    'Use only one of "symbols" and "numbers".')
            else:
                self.new_array('numbers', symbols2numbers(symbols), int)

        if cell is None:
            cell = np.eye(3)
        self.set_cell(cell)

        if positions is None:
            if scaled_positions is None:
                positions = np.zeros((len(self.arrays['numbers']), 3))
            else:
                positions = np.dot(scaled_positions, self._cell)
        else:
            if scaled_positions is not None:
                raise RuntimeError('Both scaled and cartesian positions set!')
        self.new_array('positions', positions, float, (3,))

        if pbc is None:
            pbc = False
        self.set_pbc(pbc)

    def set_cell(self, cell, scale_atoms=False, fix=None):
        if fix is not None:
            raise TypeError('Please use scale_atoms=%s' % (not fix))

        cell = np.array(cell, float)
        if cell.shape == (3,):
            cell = np.diag(cell)
        elif cell.shape != (3, 3):
            raise ValueError('Cell must be length 3 sequence or '
                             '3x3 matrix!')
        if scale_atoms:
            M = np.linalg.solve(self._cell, cell)
            self.arrays['positions'][:] = np.dot(self.arrays['positions'], M)
        self._cell = cell

    def get_cell(self):
        return self._cell.copy()

    def get_reciprocal_cell(self):
        rec_unit_cell = np.linalg.inv(self.get_cell()).transpose()
        return rec_unit_cell

    def set_pbc(self, pbc):
        if isinstance(pbc, int):
            pbc = (pbc,) * 3
        self._pbc = np.array(pbc, bool)

    def get_pbc(self):
        return self._pbc.copy()

    def new_array(self, name, a, dtype=None, shape=None):
        if dtype is not None:
            a = np.array(a, dtype)
        else:
            a = a.copy()


        if name in self.arrays:
            raise RuntimeError

        for b in self.arrays.values():
            if len(a) != len(b):
                raise ValueError('Array has wrong length: %d != %d.' %
                                 (len(a), len(b)))
            break

        if shape is not None and a.shape[1:] != shape:
            raise ValueError('Array has wrong shape %s != %s.' %
                             (a.shape, (a.shape[0:1] + shape)))

        self.arrays[name] = a

    def get_array(self, name, copy=True):
        if copy:
            return self.arrays[name].copy()
        else:
            return self.arrays[name]

    def set_array(self, name, a, dtype=None, shape=None):

        b = self.arrays.get(name)
        if b is None:
            if a is not None:
                self.new_array(name, a, dtype, shape)
        else:
            if a is None:
                del self.arrays[name]
            else:
                a = np.asarray(a)
                if a.shape != b.shape:
                    raise ValueError('Array has wrong shape %s != %s.' %
                                     (a.shape, b.shape))
                b[:] = a

    def has(self, name):
        return name in self.arrays

    def set_atomic_numbers(self, numbers):
        self.set_array('numbers', numbers, int, ())

    def get_atomic_numbers(self):
        return self.arrays['numbers'].copy()

    def set_chemical_symbols(self, symbols):
        self.set_array('numbers', symbols2numbers(symbols), int, ())

    def get_chemical_symbols(self, reduce=False):
        if reduce:
            return self.get_chemical_formula(mode='reduce')
        return [ChemicalSymbols[Z] for Z in self.arrays['numbers']]

    def get_chemical_formula(self, mode='hill'):
        if len(self) == 0:
            return ''

        if mode == 'reduce':
            numbers = self.get_atomic_numbers()
            n = len(numbers)
            changes = np.concatenate(([0], np.arange(1, n)[numbers[1:] !=
                                                           numbers[:-1]]))
            symbols = [ChemicalSymbols[e] for e in numbers[changes]]
            counts = np.append(changes[1:], n) - changes
        elif mode == 'hill':
            numbers = self.get_atomic_numbers()
            elements = np.unique(numbers)
            symbols = np.array([ChemicalSymbols[e] for e in elements])
            counts = np.array([(numbers == e).sum() for e in elements])

            ind = symbols.argsort()
            symbols = symbols[ind]
            counts = counts[ind]

            if 'H' in symbols:
                i = np.arange(len(symbols))[symbols == 'H']
                symbols = np.insert(np.delete(symbols, i), 0, symbols[i])
                counts = np.insert(np.delete(counts, i), 0, counts[i])
            if 'C' in symbols:
                i = np.arange(len(symbols))[symbols == 'C']
                symbols = np.insert(np.delete(symbols, i), 0, symbols[i])
                counts = np.insert(np.delete(counts, i), 0, counts[i])
        elif mode == 'all':
            numbers = self.get_atomic_numbers()
            symbols = [ChemicalSymbols[n] for n in numbers]
            counts = [1] * len(numbers)
        else:
            raise ValueError("Use mode = 'all', 'reduce' or 'hill'.")

        formula = ''
        for s, c in zip(symbols, counts):
            formula += s
            if c > 1:
                formula += str(c)
        return formula

    def set_positions(self, newpositions):
        positions = self.arrays['positions']

        self.set_array('positions', newpositions, shape=(3,))

    def get_positions(self, wrap=False):
        if wrap:
            scaled = self.get_scaled_positions()
            return np.dot(scaled, self._cell)
        else:
            return self.arrays['positions'].copy()

    def copy(self):
        import copy
        atoms = self.__class__(cell=self._cell, pbc=self._pbc)

        atoms.arrays = {}
        for name, a in self.arrays.items():
            atoms.arrays[name] = a.copy()

        return atoms

    def __len__(self):
        return len(self.arrays['positions'])

    def get_number_of_atoms(self):
        return len(self)

    def __repr__(self):
        num = self.get_atomic_numbers()
        N = len(num)
        if N == 0:
            symbols = ''
        elif N <= 60:
            symbols = self.get_chemical_formula('reduce')
        else:
            symbols = self.get_chemical_formula('hill')
        s = "%s(symbols='%s', " % (self.__class__.__name__, symbols)
        for name in self.arrays:
            if name == 'numbers':
                continue
            s += '%s=..., ' % name
        if (self._cell - np.diag(self._cell.diagonal())).any():
            s += 'cell=%s, ' % self._cell.tolist()
        else:
            s += 'cell=%s, ' % self._cell.diagonal().tolist()
        s += 'pbc=%s, ' % self._pbc.tolist()
        return s[:-2] + ')'

    def __add__(self, other):
        atoms = self.copy()
        atoms += other
        return atoms

#    def __getitem__(self, i):
#        if isinstance(i, int):
#            natoms = len(self)
#            if i < -natoms or i >= natoms:
#                raise IndexError('Index out of range.')
#
#            return Atom(atoms=self, index=i)

#        atoms = self.__class__(cell=self._cell, pbc=self._pbc, info=self.info)
#        atoms.adsorbate_info = self.adsorbate_info
#
#        atoms.arrays = {}
#        for name, a in self.arrays.items():
#            atoms.arrays[name] = a[i].copy()
#
#        return atoms

    def __delitem__(self, i):

        mask = np.ones(len(self), bool)
        mask[i] = False
        for name, a in self.arrays.items():
            self.arrays[name] = a[mask]

    def extend(self, other):
        if isinstance(other, Atom):
            other = self.__class__([other])

        n1 = len(self)
        n2 = len(other)

        for name, a1 in self.arrays.items():
            a = np.zeros((n1 + n2,) + a1.shape[1:], a1.dtype)
            a[:n1] = a1
            if name == 'masses':
                a2 = other.get_masses()
            else:
                a2 = other.arrays.get(name)
            if a2 is not None:
                a[n1:] = a2
            self.arrays[name] = a

        for name, a2 in other.arrays.items():
            if name in self.arrays:
                continue
            a = np.empty((n1 + n2,) + a2.shape[1:], a2.dtype)
            a[n1:] = a2
            if name == 'masses':
                a[:n1] = self.get_masses()[:n1]
            else:
                a[:n1] = 0

            self.set_array(name, a)

        return self

    __iadd__ = extend

    def append(self, atom):
        self.extend(self.__class__([atom]))

    def pop(self, i=-1):
        atom = self[i]
        atom.cut_reference_to_atoms()
        del self[i]
        return atom

    def __imul__(self, m):
        if isinstance(m, int):
            m = (m, m, m)

        M = np.product(m)
        n = len(self)

        for name, a in self.arrays.items():
            self.arrays[name] = np.tile(a, (M,) + (1,) * (len(a.shape) - 1))

        positions = self.arrays['positions']
        i0 = 0
        for m0 in range(m[0]):
            for m1 in range(m[1]):
                for m2 in range(m[2]):
                    i1 = i0 + n
                    positions[i0:i1] += np.dot((m0, m1, m2), self._cell)
                    i0 = i1

        self._cell = np.array([m[c] * self._cell[c] for c in range(3)])

        return self

    def repeat(self, rep):
        atoms = self.copy()
        atoms *= rep
        return atoms

    __mul__ = repeat

    def translate(self, displacement):
        self.arrays['positions'] += np.array(displacement)

    def get_distance(self, a0, a1, mic=False):
        R = self.arrays['positions']
        D = R[a1] - R[a0]
        if mic:
            Dr = np.linalg.solve(self._cell.T, D)
            D = np.dot(Dr - np.round(Dr) * self._pbc, self._cell)
        return np.linalg.norm(D)

    def set_distance(self, a0, a1, distance, fix=0.5):
        R = self.arrays['positions']
        D = R[a1] - R[a0]
        x = 1.0 - distance / np.linalg.norm(D)
        R[a0] += (x * fix) * D
        R[a1] -= (x * (1.0 - fix)) * D

    def get_scaled_positions(self):
        scaled = np.linalg.solve(self._cell.T, self.arrays['positions'].T).T
        for i in range(3):
            if self._pbc[i]:
                scaled[:, i] %= 1.0
                scaled[:, i] %= 1.0
        return scaled

    def set_scaled_positions(self, scaled):
        self.arrays['positions'][:] = np.dot(scaled, self._cell)

    def __eq__(self, other):
        try:
            a = self.arrays
            b = other.arrays
            return (len(self) == len(other) and
                    (a['positions'] == b['positions']).all() and
                    (a['numbers'] == b['numbers']).all() and
                    (self._cell == other.cell).all() and
                    (self._pbc == other.pbc).all())
        except AttributeError:
            return NotImplemented

    def __ne__(self, other):
        eq = self.__eq__(other)
        if eq is NotImplemented:
            return eq
        else:
            return not eq

    __hash__ = None

    def get_volume(self):
        return abs(np.linalg.det(self._cell))

    def _get_positions(self):
        return self.arrays['positions']

    def _set_positions(self, pos):
        self.arrays['positions'][:] = pos

    positions = property(_get_positions, _set_positions,
                         doc='Attribute for direct ' +
                         'manipulation of the positions.')

    def _get_atomic_numbers(self):
        return self.arrays['numbers']

    numbers = property(_get_atomic_numbers, set_atomic_numbers,
                       doc='Attribute for direct ' +
                       'manipulation of the atomic numbers.')

    def _get_cell(self):
        return self._cell

    cell = property(_get_cell, set_cell, doc='Attribute for direct ' +
                       'manipulation of the unit cell.')

    def _get_pbc(self):
        return self._pbc

    pbc = property(_get_pbc, set_pbc,
                   doc='Attribute for direct manipulation ' +
                   'of the periodic boundary condition flags.')

    def get_name(self):
        return self.get_chemical_formula(mode='hill')

def string2symbols(s):
    n = len(s)

    if n == 0:
        return []

    c = s[0]

    if c.isdigit():
        i = 1
        while i < n and s[i].isdigit():
            i += 1
        return int(s[:i]) * string2symbols(s[i:])

    if c == '(':
        p = 0
        for i, c in enumerate(s):
            if c == '(':
                p += 1
            elif c == ')':
                p -= 1
                if p == 0:
                    break
        j = i + 1
        while j < n and s[j].isdigit():
            j += 1
        if j > i + 1:
            m = int(s[i + 1:j])
        else:
            m = 1
        return m * string2symbols(s[1:i]) + string2symbols(s[j:])

    if c.isupper():
        i = 1
        if 1 < n and s[1].islower():
            i += 1
        j = i
        while j < n and s[j].isdigit():
            j += 1
        if j > i:
            m = int(s[i:j])
        else:
            m = 1
        return m * [s[:i]] + string2symbols(s[j:])
    else:
        raise ValueError

def symbols2numbers(symbols):
    if isinstance(symbols, str):
        symbols = string2symbols(symbols)
    numbers = []
    for s in symbols:
        if isinstance(s, str):
            numbers.append(AtomicNumbers[s])
        else:
            numbers.append(s)
    return numbers

def string2vector(v):
    if isinstance(v, str):
        if v[0] == '-':
            return -string2vector(v[1:])
        w = np.zeros(3)
        w['xyz'.index(v)] = 1.0
        return w
    return np.array(v, float)
