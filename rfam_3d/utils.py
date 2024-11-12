# -*- coding: utf-8 -*-

# Copyright [2009-2024] EMBL-European Bioinformatics Institute
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""A series of utilities that are used in the package.
"""

from typing import NoReturn


def asupper(s: str) -> str:
    """Upper a string.

    >>> asupper("bob")
    "BOB"
    >>> asupper("BOB")
    "BOB"
    """
    return s.upper()


def aslower(s: str) -> str:
    """Lowercase a string.

    >>> aslower("BOB")
    "bob"
    >>> aslower("bob")
    "bob"
    >>> aslower("Bob")
    "bob"
    """
    return s.lower()


def to_camel_case(snake_str: str) -> str:
    """Covert a snake case string to camel case.

    >>> to_camel_case("a_string")
    "aString"
    >>> to_camel_case("astring")
    "astring"
    """
    components = snake_str.split("_")
    return components[0] + "".join(x.title() for x in components[1:])


def assert_never(x: NoReturn) -> NoReturn:
    """A function, which is meant to serve as a type hint and a runtime
    check. Basically this is used at the end of match statements as the
    fallback branch and will cause a runtime failure if that branch is taken.
    """
    assert False, "Unhandled type: {}".format(type(x).__name__)
