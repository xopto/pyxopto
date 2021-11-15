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

import sys
import os
import re
import subprocess
import platform


MSVCVER =  [100,  110,  120,  130,  140,  150,  160 ]
MSVS =     [2010, 2012, 2013, 2015, 2015, 2017, 2019]
MS_VS2VER = {2010:100, 2012:110, 2013:120, 2015:130, 2017:150, 2019:160}
MS_VER2VS = {100:2010, 110:2012, 120:2013, 130:2015, 150:2017, 160:2019}

def os_bits() -> int:
    '''
    Returns 32 for 32-bit OS, 64 for 64-bit OS.
    '''
    if sys.maxsize > 4294967295:
        return 64
    else:
        return 32

def subprocess_cmd(command: list, verbose: bool = False, **kwargs) -> str:
    '''
    Execute a command using :py:class:`subprocess.Popen`.

    Parameters
    ----------
    command: list
        List of str that represent the command and all the arguments.
    verbose: bool
        Turn on verbose reporting.
    kwargs: dict
        Keyword arguments passed to :py:meth:`subprocess.Popen.__init__`.

    Returns
    -------
    err: str
        Error from stderr if any.
    '''
    if verbose:
        print('Executing shell command:', command)
    proc = subprocess.Popen(
        command,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, **kwargs
    )
    out, err = proc.communicate()
    if verbose:
        print('Command stdout:', out.decode('utf8'), sep='\n\t')
        print('Command stderr:', err.decode('utf8'), sep='\n\t')

    return err.decode('utf8')

def find_msvcver(supported: list = None, allver: bool = False) -> str:
    '''
    Find a version of installed MSC build tools.

    Parameters
    ----------
    allver: bool
        If True, all the found versions are returned, else only the highest
        version is returned.
    supported: list
        Optional list of supported MSC build tool versions.

    Returns
    -------
    ver: str
        A list of found MSC build tools if allver is True, else only the
        highest found version. Returns [] if no MSC build tools are found.
    '''
    if supported is None:
        supported = MSVCVER
    found = []

    t = re.compile('VS(?P<version>[0-9][0-9][0-9])COMNTOOLS')
    for var in os.environ:
        res = t.match(var)
        if res is not None:
            version = int(res.group('version'))
            if version in supported:
                found.append(version)
    if found:
        if not allver:
            found = max(found)
        return found

def msvcver2vs(msvcver):
    if msvcver in MS_VER2VS:
        return MS_VER2VS.get(msvcver)
    return 'unknown'

def find_msvsver(supported=None, allver=False):
    '''
    Find a version of installed Visual studio tools.

    Parameters
    ----------
    allver: bool
        If True, all the found versions are returned, else only the highest
        version is returned.
    supported: list
        Optional list of supported MSC build tool versions.

    Returns
    -------
    ver: str
        A list of found Visual Studio installations if allver is True, else
        only the highest found version.
        Returns [] if no MSC build tools are found.
    '''
    vers = find_msvcver(supported, allver)
    if isinstance(vers, int):
        return msvcver2vs(vers)
    else:
        msvs = []
        for ver in vers:
            msvs.append(msvcver2vs(ver))
        return msvs

def msvc_build_dll(srclist, msvcver=None, target=None, ccopts='',
                   inclist=None, deflist=None, libdirlist=None, liblist=None,
                   verbose=False, envvars=None, cwd=None, moveto=None,
                   implib=False, automachine=True, keepscript=False):
    '''
    Builds a DLL using MSC build tools.
    '''

    machineopts = ''
    if not automachine:
        machine = platform.machine().lower()
        if machine in ['x86', 'i386', 'ia-32', 'ia_32']:
            machineopts = '/MACHINE:X86'
        elif machine in ['x64', 'amd64', 'x86-64', 'x86_64',
                         'ia-32e', 'ia_32e', 'em64t', 'intel64']:
            if os_bits() > 32:
                machineopts = '/MACHINE:X64'
            else:
                machineopts = '/MACHINE:x86'
        elif machine in ['arm']:
            machineopts = '/MACHINE:ARM'
        elif machine in ['ebc']:
            machineopts = '/MACHINE:EBC'
            raise ValueError('Target platform {} not supported '\
                'by the compiller!'.format(machine))

    if verbose:
        if os_bits() > 32:
            print('Building 64-bit dll on {} machine!'.format(
                platform.machine()))
        else:
            print('Building 32-bit dll on {} machine!'.format(
                platform.machine()))

    if msvcver is None:
        msvcver = find_msvcver()
        if msvcver is None:
            raise RuntimeError(
                'Supported version of Microsoft Visual Studio tools not fund!')
    else:
        if msvcver not in find_msvcver(allver=True):
            raise RuntimeError('The specified version of '\
                'Microsoft Visual Studio {} and '\
                '"VS{}COMNTOOLS tools" not fund!'.format(
                    MS_VER2VS[msvcver], msvcver))

    if verbose:
        print('Found Microsoft Visual Studio {}... VS{}COMNTOOLS!'.format(
            msvcver2vs(msvcver), msvcver))

    path = os.environ['VS{}COMNTOOLS'.format(msvcver)]
    if os_bits() > 32:
        init_cmd = '"{}\\..\\..\\VC\\bin\\x86_amd64\\vcvarsx86_amd64.bat"'.format(path)
    else:
        init_cmd = '"{}\\..\\..\\VC\\bin\\vcvars32.bat"'.format(path)

    if target is None:
        name = os.path.basename(srclist[0])
        target = os.path.splitext(name)[0] + '.dll'


    implibname = os.path.splitext(os.path.basename(target))[0] + '.lib'

    env = None
    if envvars is not None:
        env = os.environ
        for var in envvars:
            env[var] = envvars[var]

    allIncludes = ''
    if inclist is not None:
        allIncludes = ('/I "{}" '*len(inclist)).format(*inclist)

    allLibdirs = ''
    if libdirlist is not None:
        allLibdirs = ('/LIBPATH:"{}" '*len(libdirlist)).format(*libdirlist)

    allSrc = ''
    allSrc = ('"{}" '*len(srclist)).format(*srclist)

    allObj = ''
    objlist = []
    for src in srclist:
        name = os.path.basename(src)
        objlist.append(os.path.splitext(name)[0] + '.obj')
    allObj = ' '.join(objlist)

    allDefs = ''
    if deflist is not None:
        allDefs += ('/D{} '*len(deflist)).format(*deflist)

    targetstr = '/OUT:"{}" '.format(target)

    implibstr = '/IMPLIB:"{}" '.format(implibname)

    allLibs = ''
    if liblist is not None:
        allLibs = ('"{}" '*len(liblist)).format(*liblist)

    if os_bits() > 32:
        build_cmd = 'cl /GS /GL /W3 /Gy /Zi /Gm- /O2 /sdl /fp:precise /WX- '\
            '/Zc:forScope /Gd /Oi /MD /EHsc /nologo '\
            '{} {} {} /LD /link {} /NOLOGO /DLL {} '\
            '/OPT:REF /OPT:ICF /INCREMENTAL:NO /LTCG /NXCOMPAT '\
            '/DYNAMICBASE /TLBID:1 {} {} {}'.format(
                allIncludes, allDefs, allSrc, allLibdirs,
                machineopts, implibstr, targetstr, allLibs)
    else:
        build_cmd = 'cl /GS /GL /analyze- /W3 /Gy /Zi /Gm- /O2 /sdl /fp:precise '\
        '/WX- /Zc:forScope /Gd /Oy- /Oi /MD /EHsc /nologo {} {} {} '\
        '/LD /link {} /LTCG /NXCOMPAT /DLL {} /OPT:REF '\
        '/SAFESEH /INCREMENTAL:NO /OPT:ICF /NOLOGO /TLBID:1 {} {} {}'.format(
            allIncludes, allDefs, allSrc, allLibdirs,
            machineopts, implibstr, targetstr, allLibs)

    del_obj_cmd = ''
    for obj in objlist:
        del_obj_cmd += 'del /Q {}{}'.format(obj, os.linesep)

    move_target_cmd = ''
    move_implib_cmd = ''
    move_mkdir_cmd = ''
    if moveto is not None:
        move_target_cmd = 'move "{}" "{}"'.format(
            target, os.path.join(moveto, target))
        if implib:
            move_implib_cmd = 'move "{}" "{}"'.format(
                implibname, os.path.join(moveto, implibname))
        move_mkdir_cmd = 'if not exist "{}" mkdir "{}"'.format(moveto, moveto)

    del_implib_cmd = ''
    if not implib:
        del_implib_cmd = 'del /Q {}'.format(implibname)

    bat_file = ['@echo off', "call " + init_cmd, build_cmd,
                move_mkdir_cmd, move_target_cmd, move_implib_cmd,
                del_implib_cmd, del_obj_cmd,
                'del /Q *.exp', 'del /Q *.pdb']

    bat_file = os.linesep.join(bat_file)

    bat_filename = 'build.bat'
    if cwd is not None:
        bat_filename = os.path.abspath(os.path.join(cwd, bat_filename))
    with open(bat_filename, 'w') as fid:
        fid.write(bat_file)

    if verbose:
        if cwd is not None:
            print('Preparing a temporary batch file build.bat in {}:'.format(
                cwd))
        else:
            print('Preparing a temporary batch file build.bat in {}:'.format(
                os.getcwd()))
        print('#'*80, bat_file, '#'*80, sep='\n')

    err = subprocess_cmd([bat_filename],
                         env=env, cwd=cwd, verbose=verbose)
    if err:
        raise RuntimeError('Build process failed:\n' + err)

    if not keepscript:
        err = subprocess_cmd(['del', '/Q', bat_filename],
                             env=env, cwd=cwd, shell=True, verbose=verbose)
        if err:
            raise RuntimeError('Could not delete the temporary files:\n' + err)

def _gcc_version(gcc='gcc'):
    p = subprocess.Popen([gcc, '--version'], stdout=subprocess.PIPE)
    out, err = p.communicate()
    return out.decode('ascii').lower()

def find_mingw_64(gcc='gcc') -> str:
    '''
    Returns the version of installed 32-bit MINGW compiller or None.
    '''
    p = subprocess.Popen([gcc, '-dumpmachine'], stdout=subprocess.PIPE)
    out, err = p.communicate()
    outstr = out.decode('ascii').lower()
    if ('x86_64' in outstr or 'x86-64' in outstr or 'x64' in outstr or
            'amd64' in outstr) and 'mingw' in outstr:
        return _gcc_version(), outstr
    return None

def find_mingw() -> str:
    '''
    Returns the version of installed MINGW compiller or None.
    If the underlaying OS is 32-bit, only 32-bit compilers are considered,
    else only 64-bit compilers are considered.
    '''
    if os_bits() > 32:
        return find_mingw_64()
    else:
        return find_mingw_32()

def find_mingw_32(gcc='gcc') -> str:
    '''
    Returns the version of installed 32-bit compiller or None.
    '''
    p = subprocess.Popen([gcc, '-dumpmachine'], stdout=subprocess.PIPE)
    out, err = p.communicate()
    outstr = out.decode('ascii').lower()
    if ('i686' in outstr or 'i386' in outstr or 'ia-32' in outstr or
            'ia_32' in outstr) and 'mingw' in outstr:
        return _gcc_version(), outstr
    return None

def mingw_build_dll(srclist, inclist=None, liblist=None, libdirlist=None, 
                    target=None, deflist=None, ccopts='',
                    gcc='gcc', opt='-O2', std=None, implib=False,
                    moveto=None, cwd=None, verbose=False, envvars=None,
                    keepscript=False):
    '''
    Builds a DLL using MIGW tools.
    '''
    if os_bits() > 32:
        cc = find_mingw_64()
        if cc is None:
            raise RuntimeError('Mingw gcc compiler not found!')
    else:
        cc = find_mingw_32()
        if cc is None:
            raise RuntimeError('Mingw gcc compiler not found!')

    if verbose:
        print('Found mingw gcc: {}\nMachine: {}!\n'.format(
            cc[0].strip(os.linesep), cc[1].strip(os.linesep)))

    env = None
    if envvars is not None:
        env = os.environ
        for var in envvars:
            env[var] = envvars[var]

    stdstr = ''
    if std is not None:
        stdstr = '-std=' + std

    allIncludes = ''
    allLibs = ''
    allLibdirs = ''
    allDefs = ''
    if inclist is not None:
        allIncludes = ('-I"{}" '*len(inclist)).format(*inclist)
    if liblist is not None:
        allLibs = ('-l{} '*len(liblist)).format(*liblist)
    if libdirlist is not None:
        allLibdirs = ('-L"{}" '*len(libdirlist)).format(*libdirlist)
    if deflist is not None:
        allDefs = ('-D{} '*len(deflist)).format(*deflist)

    allSrc = ('"{}" '*len(srclist)).format(*srclist)
    allObj = ''
    objlist = []
    for src in srclist:
        name = os.path.basename(src)
        objlist.append(os.path.splitext(name)[0] + '.o')
    allObj = ' '.join(objlist)

    if target is None:
        target = os.path.basename(srclist[0])
        target = os.path.splitext(target)[0] + '.dll'

    targetstr = '-o ' + target

    implibstr = ''
    implibname = ''
    if implib:
        implibname = os.path.splitext(os.path.basename(target))[0] + '.lib'
        implibstr = '-Wl,--out-implib={}'.format(implibname)

    gcc_build_cmd = 'gcc -c {} {} {} {} {} {}'.format(
        opt, ccopts, stdstr, allDefs, allIncludes, allSrc)
    gcc_link_cmd = 'gcc -shared -fPIC {} {} {} {} {}'.format(
        allLibdirs, targetstr, allObj, allLibs, implibstr)
    strip_cmd = 'strip --strip-all --discard-all {}'.format(target)

    del_obj_cmd = ''
    for obj in objlist:
        del_obj_cmd += 'del /Q {}{}'.format(obj, os.linesep)

    move_target_cmd = ''
    move_implib_cmd = ''
    move_mkdir_cmd = ''
    if moveto is not None:
        move_target_cmd = 'move "{}" "{}"'.format(
            target, os.path.join(moveto, target))
        if implib:
            move_implib_cmd = 'move "{}" "{}"'.format(
                implibname, os.path.join(moveto, implibname))
        move_mkdir_cmd = 'if not exist "{}" mkdir "{}"'.format(moveto, moveto)

    bat_file = os.linesep.join(
        ['@echo off', gcc_build_cmd, gcc_link_cmd, strip_cmd,
         move_mkdir_cmd, move_target_cmd, move_implib_cmd, del_obj_cmd])

    bat_filename = 'build.bat'
    if cwd is not None:
        bat_filename = os.path.abspath(os.path.join(cwd, bat_filename))
    with open(bat_filename, 'w') as fid:
        fid.write(bat_file)

    if verbose:
        if cwd is not None:
            print('Preparing a temporary batch file build.bat in {}:'.format(
                cwd))
        else:
            print('Preparing a temporary batch file build.bat in {}:'.format(
                os.getcwd()))
        print('#'*80, bat_file, '#'*80, sep='\n')

    err = subprocess_cmd([bat_filename], env=env, cwd=cwd, verbose=verbose)
    if err:
        raise RuntimeError('Build process failed:\n' + err)

    if not keepscript:
        err = subprocess_cmd(['del', '/Q', bat_filename],
                             env=env, cwd=cwd, shell=True, verbose=verbose)
        if err:
            raise RuntimeError('Could not delete the temporary files:\n' + err)


def find_gcc_64(gcc='gcc') -> str:
    '''
    Returns version of the installed 64-bit GCC compiler or None.
    '''
    p = subprocess.Popen([gcc, '-dumpmachine'], stdout=subprocess.PIPE)
    out, err = p.communicate()
    outstr = out.decode('ascii').lower()
    if ('x86_64' in outstr or 'x86-64' in outstr or 'x64' in outstr or
            'amd64' in outstr) and 'linux' in outstr:
        return _gcc_version(), outstr
    return None

def find_gcc_32(gcc='gcc'):
    '''
    Returns version of the installed 32-bit GCC compiler or None.
    '''
    p = subprocess.Popen([gcc, '-dumpmachine'], stdout=subprocess.PIPE)
    out, err = p.communicate()
    outstr = out.decode('ascii').lower()
    if ('i686' in outstr or 'i386' in outstr or 'ia-32' in outstr or
            'ia_32' in outstr) and 'linux' in outstr:
        return _gcc_version(), outstr
    return None

def find_gcc():
    '''
    Find a version of installed GCC that matches the number OS bits (32 or 64). 
    '''
    if os_bits() > 32:
        return find_gcc_64()
    else:
        return find_gcc_32()

def gcc_build_so_dll(srclist, inclist=None, liblist=None, libdirlist=None,
                  target=None, deflist=None,
                  gcc='gcc', opt='-O2', soname=None, std=None, implib=False,
                  moveto=None, cwd=None, verbose=False, envvars=None,
                  keepscript=False, shell='bash'):
    '''
    Builds a DLL (SO) using GCC.
    '''
    implib = False  # ... no import libraries on linux platform
    if os_bits() > 32:
        cc = find_gcc_64()
        if cc is None:
            raise RuntimeError('GNU gcc compiler not found!')
    else:
        cc = find_gcc_32()
        if cc is None:
            raise RuntimeError('GNU gcc compiler not found!')

    if verbose:
        print('Found GNU gcc: {}\nMachine: {}!\n'.format(
            cc[0].strip(os.linesep), cc[1].strip(os.linesep)))

    env = None
    if envvars is not None:
        env = os.environ
        for var in envvars:
            env[var] = envvars[var]

    stdstr = ''
    if std is not None:
        stdstr = '-std=' + std

    allIncludes = ''
    allLibs = ''
    allLibdirs = ''
    allDefs = ''
    if inclist is not None:
        allIncludes = ('-I"{}" '*len(inclist)).format(*inclist)
    if liblist is not None:
        allLibs = ('-l{} '*len(liblist)).format(*liblist)
    if libdirlist is not None:
        allLibdirs = ('-L"{}" '*len(libdirlist)).format(*libdirlist)
    if deflist is not None:
        allDefs = ('-D{} '*len(deflist)).format(*deflist)

    allSrc = ('"{}" '*len(srclist)).format(*srclist)
    allObj = ''
    objlist = []
    for src in srclist:
        name = os.path.basename(src)
        objlist.append(os.path.splitext(name)[0] + '.o')
    allObj = ' '.join(objlist)

    if target is None:
        target = os.path.basename(srclist[0])
        target = os.path.splitext(target)[0] + '.so'

    targetstr = '-o ' + target

    implibstr = ''
    implibname = ''
    if implib:
        implibname = os.path.splitext(os.path.basename(target))[0] + '.a'
        implibstr = '-Wl,--out-implib={}'.format(implibname)

    sonameStr = ''
    if soname is not None:
        sonameStr = '-Wl,-soname={}.{}'.format(target, str(soname))

    gcc_build_cmd = 'gcc -fPIC -c {} {} {} {} {}'.format(
        opt, stdstr, allDefs, allIncludes, allSrc)
    gcc_link_cmd = 'gcc -shared -fPIC {} {} {} {} {} {}'.format(
        sonameStr, allLibdirs, targetstr, allObj, allLibs, implibstr)
    strip_cmd = 'strip --strip-all --discard-all {}'.format(target)

    del_obj_cmd = ''
    for obj in objlist:
        del_obj_cmd += 'rm -f {}{}'.format(obj, os.linesep)

    move_target_cmd = ''
    move_implib_cmd = ''
    move_mkdir_cmd = ''
    if moveto is not None:
        move_target_cmd = 'mv "{}" "{}"'.format(
            target, os.path.join(moveto, target))
        if implib:
            move_implib_cmd = 'mv "{}" "{}"'.format(
                implibname, os.path.join(moveto, implibname))
        move_mkdir_cmd = 'mkdir -p "{}"'.format(moveto)

    script_file = os.linesep.join([
        '#!/bin/{}'.format(shell), '',
        '# this will make the script abort on first error',
        'set -e', 'set -o pipefail', '',
        gcc_build_cmd, gcc_link_cmd, strip_cmd,
        move_mkdir_cmd, move_target_cmd, move_implib_cmd,
        del_obj_cmd, os.linesep])

    script_filename = 'build.sh'
    if cwd is not None:
        script_filename = os.path.abspath(os.path.join(cwd, script_filename))
    with open(script_filename, 'w') as fid:
        fid.write(script_file)

    if verbose:
        if cwd is not None:
            print('Preparing a temporary batch file build.sh in {}:'.format(
                cwd))
        else:
            print('Preparing a temporary batch file build.sh in {}:'.format(
                os.getcwd()))
        print('#'*80, script_file, '#'*80, sep='\n')

    err = subprocess_cmd([shell, script_filename],
                         env=env, cwd=cwd, verbose=verbose)
    if err:
        raise RuntimeError('Build process failed:\n' + err)

    if not keepscript:
        err = subprocess_cmd(['rm -f "{}"'.format(script_filename)],
                             env=env, cwd=cwd, shell=True, verbose=verbose)
        if err:
            raise RuntimeError('Could not delete the temporary files:\n' + err)
#%%

def build_so_dll(*args, cc: str or list = None, **kwargs) -> bool:
    '''
    Build a dll or so library using the available native compiles.

    Parameters
    ----------
    args: list
        Positional parameters passed to the DLL / SO build function.
    kwargs: dict
        Keyword parameters passed to the DLL / SO build function.

    Returns
    -------
    res: bool
        True on success, else False.
    '''
    if isinstance(cc, str):
        cc = [cc]

    if 'windows' in platform.system().lower():
        if cc is None:
            cc = ['msvc', 'mingw']
        # build
    else:
        cc = ['gcc']

    done = False
    for compiler, index in zip(cc, range(len(cc))):
        if 'msvcver' in kwargs:
            kwargs.pop('msvcver')

        compiler = str(compiler).lower()
        if 'mingw' in compiler:
            build_fun = mingw_build_dll
        elif 'gcc' in compiler:
            build_fun = gcc_build_so_dll
        elif 'msvc' in compiler or 'ms' in compiler or 'msvs' in compiler:
            for version in MSVS:
                msvcver = None
                if str(version) in compiler:
                    msvcver = MS_VS2VER.get(version)
                    kwargs['msvcver'] = msvcver
                    break
            build_fun = msvc_build_dll
        else:
            raise ValueError('Compiler {} not supported!'.format(compiler))

        try:
            build_fun(*args, **kwargs)
            done = True
            break
        except RuntimeError as err:
            if index >= len(cc) - 1:
                raise err
    return done

def find_compilers():
    osname = platform.system().lower()
    index = 1
    if 'windows' in osname:
        msvcver = find_msvcver(allver=True)
        for ver in msvcver:
            print('{}. Microsoft Visual Studio {}.\n'.format(
                index, MS_VER2VS.get(ver)))
            index += 1
        mingw = find_mingw()
        if mingw is not None:
            print('{}. MINGW gcc: {}\nMachine: {}!\n'.format(
                index, mingw[0].strip(os.linesep), mingw[1].strip(os.linesep)))
            index += 1
    elif 'linux' in osname:
        gcc = find_gcc()
        if gcc is not None:
            print('{}. GNU gcc: {}\nMachine: {}!\n'.format(
                index, gcc[0].strip(os.linesep), gcc[1].strip(os.linesep)))
        index += 1


if __name__ == '__main__':
    verbose = True
    cwd = None # '/home/miran/tmp/pylibtest'
    build_target = 'rng{}.so'.format(os_bits())
    moveto = 'bin'

    find_compilers()


    if platform.system().lower() == 'windows':
        build_so_dll(
            srclist=['src/rng.cpp'],
            deflist=['RNG_EXPORTS'],
            target='rng64.dll',
            cc='mingw',
            verbose=1
        )

        '''
        build_dll(srclist=['src/rng.cpp'], inclist=['c:/windows', 'd:/tmp'],
                    libdirlist=['c:/lib', 'd:/none'], deflist=['RNG_EXPORTS'],
                    target=build_target, moveto=moveto,
                    cwd=cwd, verbose=verbose, implib=False, keepscript=True)

        build_dll(srclist=['src/rng.cpp'], inclist=['c:/windows', 'd:/tmp'],
                    libdirlist=['c:/lib', 'd:/none'], deflist=['RNG_EXPORTS'],
                    target=build_target, moveto=moveto,
                    cwd=cwd, verbose=verbose, implib=False, keepscript=True)
        '''

    if platform.system().lower() == 'linux':
        build_so_dll(
            srclist=['src/rng.cpp'],
            deflist=['RNG_EXPORTS'],
            target='rng64.so',
            cc='mingw',
            verbose=1
        )

        gcc_build_so_dll(
            srclist=['src/rng.cpp'],
            inclist=['/home/miran/src', '/home/miran/tmp'],
            libdirlist=['/home/miran/opt', '/home/miran/src'],
            deflist=['RNG_EXPORTS'],
            target=build_target,
            moveto=moveto,
            soname='1',
            cwd=cwd,
            verbose=verbose,
            keepscript=True
        )

        gcc_build_so_dll(
            srclist=['src/rng.cpp'],
            target='rng64.so',
            moveto='bin',
            verbose=1
        )
