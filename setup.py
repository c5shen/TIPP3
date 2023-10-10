#!/usr/bin/env python

# ##########################################################################
#    Copyright 2012 Siavash Mirarab, Nam Nguyen, and Tandy Warnow.
#    This file is part of SEPP.
#
#    SEPP is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    SEPP is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with SEPP.  If not, see <http://www.gnu.org/licenses/>.
# ##########################################################################

import os
import platform
import sys

# from distutils.core import setup
from distribute_setup import use_setuptools
import shutil
from setuptools import find_packages
from distutils.core import setup, Command
from distutils.spawn import find_executable
import sepp

use_setuptools(version="0.6.24")
version = "1.1.0"


def get_tools_dir(where):
    path = os.path.join(os.getcwd(), "tools", where)
    if not os.path.exists(path):
        raise OSError("TIPP does not bundle tippJsonMerger for '%s' at this time!" )
    return path


def get_tool_name(tool, bits):
    if platform.system() == "Darwin" or not bits:  # MAC doesn't have 32/64
        return tool
    is_64bits = sys.maxsize > 2**32
    return "%s-%s" % (tool, "64" if is_64bits else "32")


class ConfigSepp(Command):
    """setuptools Command"""
    description = "Configures Sepp for the current user"
    user_options = [(
        'contained', 'c',
        ("Whether SEPP should be installed in a self-contained "
         "manner or on user's home"))]

    def initopts(self):
        self.contained = None
        self.configfile = None
        self.basepath = None

    def initpath(self, name):
        if self.contained:
            self.configfile = os.path.expanduser(
                os.path.abspath(os.path.join(".tipp", name)))
            self.basepath = os.path.dirname(self.configfile)
        else:

            self.configfile = os.path.expanduser("~/.tipp/%s" % name)
            self.basepath = os.path.expanduser("~/.tipp")
        with open('home.path', 'w') as fo:
            fo.write(self.basepath)
            fo.close()

    def get_tools_dest(self):
        return os.path.join(self.basepath, "bundled-v%s" % version)

    def copy_tool_to_lib(self, tool, where, bits=True):
        shutil.copy2(
            os.path.join(get_tools_dir(where), get_tool_name(tool, bits)),
            os.path.join(self.get_tools_dest(), tool))

    def initialize_options(self):
        """init options"""
        self.initopts()

    def finalize_options(self):
        """finalize options"""
        self.initpath("main.config")
        print("\nCreating main sepp config file at %s and tools at %s" % (
            self.configfile, self.basepath))

    def run(self):
        def get_tool_name(tool, bits):
            # MAC doesn't have 32/64
            if platform.system() == "Darwin" or not bits:
                return tool
            is_64bits = sys.maxsize > 2**32
            return "%s-%s" % (tool, "64" if is_64bits else "32")

        # Create the default config file
        if not os.path.exists(self.basepath):
            os.mkdir(self.basepath)
        if not os.path.exists(self.get_tools_dest()):
            os.mkdir(self.get_tools_dest())
        c = open("default.main.config")
        d = open(self.configfile, "w")
        for l1 in c:
            l1 = l1.replace("~", self.get_tools_dest())
            d.write(l1)
        d.close()

        # Copy tools to a bundled directory inside .tipp
        self.copy_tool_to_lib("guppy")
        self.copy_tool_to_lib("pplacer")
        self.copy_tool_to_lib("hmmalign")
        self.copy_tool_to_lib("hmmsearch")
        self.copy_tool_to_lib("hmmbuild")
        # TODO: should we compile and build merge.jar?
        self.copy_tool_to_lib("seppJsonMerger.jar", where="merge", bits=False)

class ConfigTIPP(ConfigSepp):
    """setuptools Command"""
    description = "Configures TIPP for the current user"
    user_options = [('contained', 'c',
                     ("Whether TIPP should be installed in a self-contained "
                      "manner or on user's home"))]

    def initialize_options(self):
        """init options"""
        self.initopts()

    def finalize_options(self):
        """finalize options"""
        self.initpath("tipp.config")
        print("\nCreating main TIPP config file at %s and tools at %s" % (
            self.configfile, self.basepath))

    def run(self):
        ## We are going to read sepp config file and write it here 
        root_p = open(os.path.join(os.path.split(os.path.split(sepp.__file__)[0])[0], "home.path")).readlines()[0].strip()
        sepp_config_path = os.path.join(root_p, "main.config")

        # Create the default config file
        if not os.path.exists(self.basepath):
          os.mkdir(self.basepath)
        if not os.path.exists(self.get_tools_dest()):
          os.mkdir(self.get_tools_dest())
        c = open(sepp_config_path)
        d = open(self.configfile, "w")
        for l3 in c:
            #This is not needed as we are reading sepp config and not default config distributed with 
            # l3 = l3.replace("~", self.get_tools_dest())
            if (l3.find('seppJsonMerger.jar') != -1):
              l3 = "path="+self.get_tools_dest()+'/tippJsonMerger.jar\n'
            d.write(l3)
        if not os.getenv('SATE') is None:
            d.write('\n[sate]\npath=%s' % os.getenv('SATE'))
        if os.getenv('BLAST') is None and not find_executable("blastn"):
            print("\nWarning! BLAST variable is not defined.  If you plan to "
                  "run TIPP for abundance profiling,"
                  " then have BLAST pointed to blastn executable. "
                  "You can also change your config to point to"
                  " blastn by including the following line in your"
                  " config:\n[blast]\npath=/location/of/blast_directory/"
                  "blastn\n")
            d.write('\n[blast]\npath=None\n')
        else:
            blastn_path = find_executable("blastn") if os.getenv('BLAST') \
                is None else os.getenv('BLAST')
            d.write('\n[blast]\npath=%s\n' % blastn_path)

        if os.getenv('REFERENCE') is None:
            print("\nWarning! REFERENCE variable is not defined. "
                  "If you plan to run TIPP for abundance profiling, then have "
                  "REFERENCE pointed to Reference directory.  You can also "
                  "change your config to point to the Reference directory by "
                  "including the following line in your config:\n[reference]"
                  "\npath=/location/of/reference_directory/\n")
        d.write('\n[reference]\npath=%s\n' % os.getenv('REFERENCE'))
        d.write('\n[tipp]\npushdown = true\n')
        d.close()

        # Copy tools to a bundled directory inside .tipp
        self.copy_tool_to_lib("tippJsonMerger.jar", where="merge", bits=False)


setup(name="tipp",
      version=version,
      description="SATe enabled phylogenetic placement.",
      packages=find_packages(),

      url="https://github.com/smirarab/sepp",
      author="Siavash Mirarab and Nam Nguyen",
      author_email="smirarab@gmail.com, namphuon@cs.utah.edu",

      license="General Public License (GPL)",
      install_requires=["dendropy >= 4.0.0", "sepp"],
      provides=["tipp"],
      scripts=["run_abundance.py","run_tipp.py","run_tipp_tool.py"],
      cmdclass={"tipp": ConfigTIPP},
      data_files=[('', ['home.path'])],

      classifiers=["Environment :: Console",
                   "Intended Audience :: Developers",
                   "Intended Audience :: Science/Research",
                   ("License :: OSI Approved :: GNU General Public "
                    "License (GPL)"),
                   "Natural Language :: English",
                   "Operating System :: OS Independent",
                   "Programming Language :: Python",
                   "Topic :: Scientific/Engineering :: Bio-Informatics"])
