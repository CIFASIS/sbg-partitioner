#!/bin/bash

# Install pre-commit hooks.

rm -rf ../.git/hooks/pre-commit

ln -rsvf pre-commit ../.git/hooks/
