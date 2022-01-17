_ALLOWED_BUILD_OPTIONS = (
    '-cl-mad-enable',
    '-cl-no-signed-zeros',
    '-cl-unsafe-math-optimizations',
    '-cl-finite-math-only',
    '-cl-fast-relaxed-math'
)

class ClBuildOption:
    def __init__(self, value: str):
        '''
        Creates an instance of one OpenCL build option. The value of the
        options is check against the allowed list of option given by the
        OpenCL standard:

        https://www.khronos.org/registry/OpenCL/sdk/1.0/docs/man/xhtml/clBuildProgram.html

        Parameters
        ----------
        value: str
            Value of the build option. Must be one of:

              - "-cl-mad-enable",
              - "-cl-no-signed-zeros",
              - "-cl-unsafe-math-optimizations",
              - "-cl-finite-math-only",
              - "-cl-fast-relaxed-math"

        Note
        ----
        Use one of the predefined OpenCL build option instances instead of
        using this class directly:

          - :py:data:`MadEnable`
          - :py:data:`NoSignedZeros`
          - :py:data:`UnsafeMathOptimizations`
          - :py:data:`FiniteMathOnly`
          - :py:data:`FastRelaxedMath`

        '''
        if value not in _ALLOWED_BUILD_OPTIONS:
            raise ValueError(
                'OpenCL build option "{}" is not valid!\n'
                'Use one of: \n"'
                '  "-cl-mad-enable", '
                '  "-cl-no-signed-zeros", '
                '  "-cl-unsafe-math-optimizations", '
                '  "-cl-finite-math-only", or '
                '  "-cl-fast-relaxed-math".'.format(
                    value, _ALLOWED_BUILD_OPTIONS))

        self._value = str(value)

    def _get_value(self) -> str:
        return self._value
    value = property(_get_value, None, None,
                     'Value of this OpenCL build option.')

    def __str__(self):
        return self._value

    def __repr__(self):
        return 'ClBuildOption("{:s}")'.format(self._value)

MadEnable = ClBuildOption('-cl-mad-enable')
'''
Allow a * b + c to be replaced by a mad. The mad computes a * b + c with
reduced accuracy. For example, some OpenCL devices implement mad as truncate
the result of a * b before adding it to c.

Source of the above  information:
https://www.khronos.org/registry/OpenCL/sdk/1.0/docs/man/xhtml/clBuildProgram.html 
'''

NoSignedZeros = ClBuildOption('-cl-no-signed-zeros')
'''
Allow optimizations for floating-point arithmetic that ignore the signedness
of zero. IEEE 754 arithmetic specifies the behavior of distinct
+0.0 and -0.0 values, which then prohibits simplification of expressions
such as x+0.0 or 0.0*x (even with "-clfinite-math" only). This option implies
that the sign of a zero result isn't significant.

Source of the above information:
https://www.khronos.org/registry/OpenCL/sdk/1.0/docs/man/xhtml/clBuildProgram.html 
'''

UnsafeMathOptimizations = ClBuildOption('-cl-unsafe-math-optimizations')
'''
Allow optimizations for floating-point arithmetic that (a) assume that
arguments and results are valid, (b) may violate IEEE 754 standard and
(c) may violate the OpenCL numerical compliance requirements as defined
in section 7.4 for single-precision floating-point, section 9.3.9
for double-precision floating-point, and edge case behavior in section 7.5.
This option includes the "-cl-no-signed-zeros" and "-cl-mad-enable options". 

Source of the above information:
https://www.khronos.org/registry/OpenCL/sdk/1.0/docs/man/xhtml/clBuildProgram.html 
'''


FiniteMathOnly = ClBuildOption('-cl-finite-math-only')
'''
Allow optimizations for floating-point arithmetic that assume that arguments
and results are not NaNs or ±∞. This option may violate the OpenCL numerical
compliance requirements defined in in section 7.4 for single-precision
floating-point, section 9.3.9 for double-precision floating-point, and edge
case behavior in section 7.5.

Source of the above information:
https://www.khronos.org/registry/OpenCL/sdk/1.0/docs/man/xhtml/clBuildProgram.html 
'''


FastRelaxedMath = ClBuildOption('-cl-fast-relaxed-math')
'''
Sets the optimization options "-cl-finite-math-only" and
"-cl-unsafe-math-optimizations". This allows optimizations for floating-point
arithmetic that may violate the IEEE 754 standard and the OpenCL numerical
compliance requirements defined in the specification in section 7.4
for single-precision floating-point, section 9.3.9 for double-precision
floating-point, and edge case behavior in section 7.5. This option causes
the preprocessor macro __FAST_RELAXED_MATH__ to be defined in the
OpenCL program.

Source of the above information:
https://www.khronos.org/registry/OpenCL/sdk/1.0/docs/man/xhtml/clBuildProgram.html 
'''