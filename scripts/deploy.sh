#!/bin/bash

if [ "$1" == "docs" ]; then
    git config --global user.email "travis@travis-ci.org"
    git config --global user.name "Travis CI"
    git clone --branch gh-pages-test https://${GH_TOKEN}@github.com/erdc/proteus-docs.git ./docs/build > /dev/null 2>&1
    if [[ -d ./docs/build ]]
    then
      make docs
      cd ./docs/build
      git add . *
      git commit -m "Travis automatic docs build: $TRAVIS_BUILD_NUMBER"
      git push origin gh-pages-test > /dev/null 2>&1
      cd ../..
    else
      echo "Error: Travis could not clone docs repository"
      exit 1
    fi
fi
