#!/usr/bin/env
#
# Primitive linear algebra tools
#
# Created by: dadwyer@lbl.gov 2013/03/12
# Last Modified: dadwyer@lbl.gov 2013/03/14
# 
# Note: Designed so that it can be replaced by Numpy, if installed.
#

from types import *
from array import array as py_array
from copy import copy

supportedTypes = [IntType, LongType, FloatType, ComplexType]

class array(object):
    "Primitive linear algebra object (vector or matrix)"
    def __init__(self, *args):
        "Constructor"
        # Check constructor type
        if len(args) == 1:
            arg = args[0]
            if type(arg) == ListType:
                self._init_by_list(arg)
            else:
                assert False, "Invalid constructor for type '%s'" % type(arg)
        else:
            self._init_by_type_and_shape(*args)

    def _init_by_list(self, argList):
        "Construct array from python list"
        # Check that argument list is properly formatted
        [shape, dtype] = self._checkArgs(argList)
        self._init_by_type_and_shape(dtype, shape)
        # Load input data into buffer
        self._setData(argList, self._ndim, 0)
        
    def _init_by_type_and_shape(self, *args):
        "Construct new array by type and shape"
        assert len(args)>1, "Invalid number of arguments: %d" % len(args)
        dtype = args[0]
        shape = args[1]
        data = py_array('i')
        assert type(dtype) == TypeType, "First argument not a type: %s" % type(dtype)
        assert type(shape) == ListType, "Second argument not a list: %s" % type(shape)
        if len(args)>2:
            assert isinstance(args[2], py_array), "Third argument not an array: %s" % type(data)
            data = copy(args[2])
        assert all( type(i) is IntType for i in shape), "Shape must be integers only: %s" % (shape)
        # Construct shape
        self._ndim = len(shape)
        self._shape = copy(shape)
        self._dtype = dtype
        self._strides = self._makeStrides( self._shape )
        # Construct internal data buffer
        size = 1
        for dim in shape: size *= dim
        self._size = size
        self._getValue = self._getValueReal
        self._setValue = self._setValueReal
        self._bufsize = self._size
        # Handle complex arrays
        if self._dtype is ComplexType:
            self._bufsize = 2 * self._size
            self._getValue = self._getValueComplex
            self._setValue = self._setValueComplex
        # Initialize data
        if len(data)<self._bufsize:
            # Extend data array if needed
            # (Needed when none provided, or casting to complex array)
            data.extend( [0]*(self._bufsize - len(data)) )
        if self._dtype is IntType:
            self._data = py_array('i', data)
        elif self._dtype is LongType:
            self._data = py_array('l', data)
        elif self._dtype is FloatType:
            self._data = py_array('d', data)
        elif self._dtype is ComplexType:
            self._data = py_array('d', data)
        else:
            assert False, "Invalid data type: %s" % self._dtype
        return

    def _checkArgs(self, argList):
        "Check the argument list structure"
        assert type(argList) is ListType, "Initializer is not a list: %s" % argList
        argListRef = argList
        shape = []
        dtype = None
        # Determine dimension and type
        while True:
            curType = type( argListRef )
            if curType is not ListType:
                # Hit an element, check type
                assert curType in supportedTypes, "Type %s is not supported." % curType
                dtype = curType
                break
            # Check list length
            curLen = len( argListRef )
            shape.append( curLen )
            # Move to next dimension
            argListRef = argListRef[0]
        # Confirm dimension and type
        checkList( argList, shape, dtype )    
        return [shape, dtype]

    def _makeStrides(self, shape):
        "Generate strides for later quick lookup"
        strides = [1]
        for dim in reversed(shape[1:]):
            strides.append( dim * strides[-1] )
        strides.reverse()
        return strides

    def _getValueReal(self, offset):
        "Private function to get value by internal index"
        return self._data[offset]

    def _getValueComplex(self, offset):
        "Private function to get value by internal index"
        return complex(self._data[offset], self._data[offset+self._size])

    def _setValueReal(self, value, offset):
        "Private function to set value by internal index"
        self._data[offset] = value

    def _setValueComplex(self, value, offset):
        "Private function to set value by internal index"
        self._data[offset] = value.real
        self._data[offset + self._size] = value.imag
            
    def _setData(self, argList, ndim, offset):
        "Recursively initialize data, based on input values"
        if ndim == 0:
            # Hit actual element
            self._setValue(argList, offset)
        else:
            localOffset = 0
            curStride = self._strides[self._ndim - ndim]
            for element in argList:
                self._setData( element, ndim-1, offset+localOffset )
                localOffset += curStride
        return

    def _asList(self):
        "Convert to python nested list type"
        return self._asSubList(self._ndim, 0)
        
    def _asSubList(self, ndim, offset):
        "Recursively convert a subcomponent of array to a python list"
        if ndim == 0:
            return self._getValue(offset)
        lst = []
        negDim = self._ndim - ndim
        localOffset = 0
        for idx in range( self._shape[negDim] ):
            lst.append( self._asSubList(ndim-1, offset+localOffset) )
            localOffset += self._strides[negDim]
        return lst
               
    def __str__(self):
        "Convert array to string"
        return self._asSubString(self._ndim, 0, 0)

    def _asSubString(self, ndim, offset, nSpace):
        "Recursively convert a subcomponent of array to a python string"
        retStr = ''
        if ndim == 0:
            retStr += '[]'
        elif ndim == 1:
            retStr += str(self._asSubList(ndim, offset)).strip(',')
        else:
            negDim = self._ndim - ndim
            localOffset = 0
            if nSpace > 0:
                retStr += ' '*(nSpace-1)
            retStr += '['
            curShape = self._shape[negDim]
            for idx in range( self._shape[negDim] ):
                if idx > 0: 
                    retStr += ' '*(nSpace+1)
                retStr += self._asSubString(ndim-1, 
                                            offset+localOffset, 
                                            nSpace+1)
                if idx == curShape-1:
                    retStr += ']'
                else:
                    retStr += '\n'
                localOffset += self._strides[negDim]
        return retStr

    def _checkIndex(self, index):
        "Check that index is valid for this array"
        if self._ndim == 1:
            assert type(index) is IntType,"1-D array must be called using m[i] notation"
            index = (index,)
        else:
            assert type(index) is TupleType,"Multidimensional array must be called using m[i,j,k] notation"
            assert len(index)==self._ndim,"Calling %d-dim array with wrong dim: %d" % (self._ndim, len(index))
        assert all( (type(i) is IntType and index[i]>=0 and index[i]<self._shape[i]) for i in range(len(index))),"Invalid index %s for array of dimension %s" % (index, self._shape)
        return

    def _getOffset(self, index):
        "Get offset for value by array index"
        if self._ndim == 1: return index
        return sum( [index[i]*self._strides[i] for i in range(self._ndim)] )

    def __getitem__(self, key):
        "Get item from array"
        self._checkIndex(key)
        return self._getValue( self._getOffset(key) )

    def __setitem__(self, key, value):
        "Set item in array"
        self._checkIndex(key) 
        return self._setValue( value, self._getOffset(key) )

    # Array/Matrix math operations
    def __add__(self, other):
        "Add two arrays (element-wise), or add value to all elements of array"
        if type(other) in supportedTypes:
            # Add a simple type to all entries of array
            return self._addSimple(other)
        elif isinstance(other, array):
            # Add two arrays
            return self._addArray(other)
        else:
            assert False,"Attempting to add unsupported type: %s" % type(other)
        return None # Can't reach here

    def __radd__(self, other):
        "Reflected addition (A+B => B+A)."
        return self.__add__(other)

    def _addSimple(self, other):
        "Add constant value to all elements"
        newType = type(self._dtype() + other)
        retArr = array(newType, self._shape, self._data)
        otherR = 0
        otherI = 0
        if type(other) == ComplexType:
            otherR = other.real
            otherI = other.imag
        else:
            otherR = other
        # Add real component, if needed
        if otherR:
            for index in range(self._size):
                retArr._data[index] += otherR
        # Add complex component, if needed
        if otherI:
            for index in range(self._size, self._size*2):
                retArr._data[index] += otherI
        return retArr

    def _addArray(self, other):
        "Add two arrays (element-wise)"
        assert self._shape == other._shape, "Array addition requires identical dimensions: %s != %s" % (self._shape, other._shape)
        newType = type(self._dtype() + other._dtype())
        retArr = array(newType, self._shape, self._data)
        for index in range(self._size):
            retArr._data[index] += other._data[index]
        # If complex, also add complex components
        if other._dtype == ComplexType:
            for index in range(self._size, self._size*2):
                retArr._data[index] += other._data[index]
        return retArr

    def __sub__(self, other):
        "Subtract two arrays (element-wise), or subtract value from all elements of array"
        return self.__add__(-1 * other)

    def __mul__(self, other):
        "Multiply two arrays, or multiply an array by a constant value"
        if type(other) in supportedTypes:
            # Multiply all entries of array by a simple type
            return self._multiplySimple(other)
        ## Disabled to maintain compatibility with Numpy
        #elif isinstance(other, array):
        #    # Multiply two arrays
        #    return self._multiplyArray(other)       
        else:
            assert False,"Attempting to multiply unsupported type: %s" % type(other)
        return None # Can't reach here

    def __rmul__(self, other):
        "Reflected multiplication (A*B => B*A)."
        return self.__mul__(other)

    def _multiplySimple(self, other):
        "Multiply all elements by a constant value"
        newType = type(self._dtype() + other)
        retArr = array(newType, self._shape, self._data)
        if newType != ComplexType:
            # Multiply real components
            for index in range(self._size):
                retArr._data[index] *= other
        else:
            # Multiply complex elements
            otherR = 0
            otherI = 0
            if type(other) == ComplexType:
                otherR = other.real
                otherI = other.imag
            else:
                otherR = other
            for index in range(self._size):
                cIdx = self._size + index
                retArr._data[index] = self._data[index]*otherR
                retArr._data[cIdx] = self._data[index]*otherI
                if self._dtype is ComplexType:
                    retArr._data[index] -= self._data[cIdx]*otherI
                    retArr._data[cIdx] += self._data[cIdx]*otherR
        return retArr

    """  # Disabled to maintain compatibility with Numpy
    def _multiplyArray(self, other):
        "Multiply arrays, following standard matrix mathematics"
        return dot(self, other)
    """

    def conj(self):
        "Return complex conjugate of this array"
        retArr = array(self._dtype, self._shape, self._data)        
        if self._dtype is ComplexType:
            for index in range(self._size, self._size*2):
                retArr._data[index] *= -1
        return retArr

    def transpose(self):
        "Return transpose of this array"
        retArr = None 
        if self._ndim == 1:
            # 1-D transpose returns self
            retArr = array(self._dtype, self._shape, self._data)
        elif self._ndim == 2:
            newShape = self._shape[:]
            newShape.reverse()
            retArr = array(self._dtype, newShape)
            for index in range(self._size):
                nR = index / self._strides[0]
                nC = index % self._strides[0]
                transIdx = nC*retArr._strides[0] + nR
                retArr._data[transIdx] = self._data[index]
                if self._dtype == ComplexType:
                    retArr._data[self._size+transIdx] = self._data[self._size+index]
        else:
            assert False,"Transpose not yet implemented for dim > 2 arrays." 
        return retArr
        

def dot(a, b):
    "Dot-product of two multi-dimensional arrays"
    assert isinstance(a,array) and isinstance(b, array),"Dot product is only valid between two arrays, not %s and %s" % (type(a), type(b)) 
    assert a._ndim > 0 and b._ndim > 0, "Attempt to multiply 0-D matrix"
    assert a._shape[-1] == b._shape[0], "Attempt to multiply matrices of incompatible dimensions: %s and %s" % (a._shape, b._shape)
    # Dimensions and shape of new matrix
    multDim = a._shape[-1]
    newShape = a._shape[:-1] + b._shape[1:]
    if a._ndim == 1:
        # No more dimensions, dot product is a single number
        retVal = 0
        for i in range( multDim ):
            retVal += a[i]*b[i]
        return retVal
    elif all( i==1 for i in newShape ):
        # No more dimensions, dot product is a single number
        aIdxs = (0,)*(a._ndim-1)
        bIdxs = (0,)*(b._ndim-1)
        retVal = 0
        for i in range( multDim ):
            retVal += (a.__getitem__( (aIdxs + (i,)) ) 
                       * b.__getitem__( ((i,) + bIdxs)) ) 
        return retVal
    # Return a new array
    newType = type(a._dtype() + b._dtype())
    retArr = array(newType, newShape)
    newDim = retArr._ndim
    curIdx = (0,) * newDim
    isDone = False
    while True:
        # Check if dimension is complete, reset indices
        curDim = newDim - 1
        while newShape[curDim] == curIdx[curDim]:
            if curDim == 0:
                isDone = True
                break  #Done
            # Reset current index into new array
            curIdx = (curIdx[:curDim-1] 
                      + (curIdx[curDim-1]+1,)
                      + (0,)*(newDim-curDim)) 
            # Check next-higher dimension
            curDim -= 1
        if isDone: break
        #print " idx=",curIdx
        newValue = 0
        aIdx = curIdx[:a._ndim-1]
        bIdx = curIdx[-(b._ndim-1):]
        for i in range( multDim ):
            # Take dot product over this dimension
            newValue += (a.__getitem__( (aIdx + (i,)) )
                         * b.__getitem__( ((i,) + bIdx) ) )
        retArr.__setitem__(curIdx,newValue)
        # Increment to next element index
        curIdx = curIdx[:-1] + (curIdx[-1] + 1,) 
    return retArr

def checkList(argList, shape, dtype):
    """ Helper function to confirm proper argument list formatting """
    if len(shape) == 0:
        # Hit actual elements
        assert type(argList+dtype()) is dtype, "Element has invalid type: %s != type( %s )" % (dtype, argList)
    else:
        assert len(argList) == shape[0], "Incorrect length: %d != len( %s )" % (shape[0], argList) 
        for element in argList:
            checkList( element, shape[1:], dtype )
    return

