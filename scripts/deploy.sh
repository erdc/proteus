#!/bin/bash


if [ "$1" == "docs" ]; then
    if [[ -d ./docs/build ]]; then
      echo "Error: ./docs/build already exists, cannot clone docs repo there; aborting"
      exit 1
    fi
    git config --global user.email "travis@travis-ci.org"
    git config --global user.name "Travis CI"
    git clone --branch gh-pages-test https://${GH_TOKEN}@github.com/erdc/proteus-docs.git ./docs/build > /dev/null 2>&1
    # check that docs folder exists
    if [[ -d ./docs/build ]]
    then
      make docs
      cd ./docs/build
      # check that we are in right repository
      if [[ $(git remote get-url origin) == "https://${GH_TOKEN}@github.com/erdc/proteus-docs.git" ]]
      then
        git add . *
        git commit -m "Travis automatic docs build: $TRAVIS_BUILD_NUMBER"
        git push origin gh-pages-test > /dev/null 2>&1
      else
        echo "Error: wrong repository for docs; aborting"
        exit 1
      fi
      cd ../..
    else
      echo "Error: Travis could not clone docs repository; aborting"
      exit 1
    fi
fi
