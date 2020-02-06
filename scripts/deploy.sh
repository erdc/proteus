#!/bin/bash

if [ "$1" == "docs" ]; then
    make docs
    git config --global user.email "travis@travis-ci.org"
    git config --global user.name "Travis CI"
    git clone https://${GH_TOKEN}@github.com/erdc/proteus-docs.git ./proteus-docs > /dev/null 2>&1
    cd ./proteus-docs
    git checkout gh-pages-test
    cp -r ../docs/build/* .
    git add . *
    git commit -m "Travis automatic docs build: $TRAVIS_BUILD_NUMBER"
    git push origin gh-pages-test > /dev/null 2>&1
    cd ..
fi
