__author__ = 'davidmurphy'


class DictLike(object):
    """
    A base class for various data structures used in the inference.
    The insertion and retrieval of attributes has OrderedDict-like properties.
    """
    def __init__(self):
        """list of attribute names in order of initialization"""
        self._keys = []

    def __contains__(self, item):
        return item in self._keys

    def __setitem__(self, key, value):
        self.__dict__[key] = value

    def __getitem__(self, item):
        return self.__dict__[item]

    def __setattr__(self, key, value):
        self[key] = value  # write attribute values to internal dictionary
        if key not in self._keys + ['_keys']:
            self._keys.append(key)  # unique, ordered list of attribute keys

    def __iter__(self):
        for k in self._keys:
            yield k

    def __str__(self):
        """return the class type as a string"""
        return self.__class__.__name__

    def copy(self):
        """
        Create a new RunStruct with all fields set to default. Copy all member variable values from current
        struct into new struct and then return new struct. If the current struct was initialized with an outdated init
        file, the copy will include the newer member variables at their default settings.
        :return copy: a copy of the current struct
        """
        copy = self.__class__()
        for (k, v) in self.items:
            copy[k] = v
        return copy

    def setitem(self, key, value, safe=True, verbose=False):
        """
        Set attributes -- "safe" option skips uninitialized keys.

        :param key: attribute name
        :param value: attribute value
        :param safe: flags to skip unrecognized attributes
        :param verbose: flags to print a warning messaage for each
        unrecognized attribute that is added
        """
        if key in self:
            # add recognized Attribute
            setattr(self, key, value)
        else:
            if not safe:
                # add unrecognized Attribute by forcing
                setattr(self, key, value)
                if verbose:
                    # print warning about the unrecognized Attribute
                    print('WARNING: Attribute "{}" not recognized but added anyway...'.format(key))

    def setitems(self, items, safe=True):
        """
        Update each member variable in with the (key, value) pairs in items using setitem function
        :param items: (key, value) pairs of class attributes
        :param safe: a flag that indicates safe use of 'setitem' -> only existining attributes can be updated
        """
        for (k, v) in items:
            self.setitem(key=k, value=v, safe=safe)

    @property
    def keys(self):
        """return the attribute names, which are the keys to the internal dictionary"""
        return self._keys

    @property
    def values(self):
        """return the attribute values, which are the values of the internal dictionary"""
        return [self[k] for k in self]

    @property
    def items(self):
        """return the attribute name and value pairs, which are the (key, value) tuples of the internal dictionary"""
        return [(k, self[k]) for k in self]

    @property
    def dict(self):
        """return dictionary of the attribute names and values"""
        return dict(self.items)
