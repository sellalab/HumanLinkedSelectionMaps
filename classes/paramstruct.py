__author__ = 'davidmurphy'

#
# class ParamVector:
#     """a class for organizing linked selection parameters"""
#     def __init__(self, bn, cn):
#         self._params = np.zeros(shape=bn+cn+1)
#         self._bidx = np.arange(bn)
#         self._cidx = np.arange(cn) + bn
#         self._tidx = bn+cn
#
#         # fixed parameters
#         self._ufix = 7.4e-8
#         self._umax = 5e-7
#         self._min_red = 1e-5
#         self._min_pi0 = 1e-5
#         self._max_pii = 0.99
#
#     def __str__(self):
#         return ' '.join(self._params.astype(str))
#
#     def __len__(self):
#         return self._params.size
#
#     def __iter__(self):
#         for i in self._params:
#             yield i
#
#     def __setitem__(self, idx, value):
#         self._params[idx] = value
#
#     def __getitem__(self, idx):
#         return self._params[idx]
#
#     @property
#     def params(self):
#         return self._params
#
#     @params.setter
#     def params(self, new_params):
#         assert len(new_params) == len(self)
#         self._params = new_params


class Param(object):
    """a class for organizing linked selection parameters"""
    def __init__(self, p0, pmin, pmax, coef=1e-2, u0=1e-8, umax=3e-8):
        # initial param value
        self._p0 = self._pcur = p0

        # min and max values allowed for param
        self._pmin = pmin
        self._pmax = pmax

        # optional attributes: u deleterious, selection coefficient
        self._u0 = self._ucur = u0
        self._umax = umax
        self._coef = coef

    @property
    def pcur(self):
        return self._pcur

    @pcur.setter
    def pcur(self, pnew):
        self._pcur = pnew

    @property
    def ucur(self):
        return self._ucur

    @ucur.setter
    def ucur(self, new_ucur):
        self._ucur = new_ucur

    @property
    def p0(self):
        return self._p0

    @property
    def pmin(self):
        return self._pmin

    @property
    def pmax(self):
        return self._pmax

    @property
    def u0(self):
        return self._u0

    @property
    def umax(self):
        return self._umax

    @property
    def coef(self):
        return self._coef

    @property
    def uratio(self):
        return self.ucur / self.u0


