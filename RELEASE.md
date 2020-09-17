- To release a new version of proteus
Update proteus/__init__.py (remove '.dev0' and increment version)
Update setup.py (remove 'dev0' and increment version)
rebuild docs and check docs
  cd .. 
  git clone git@github.com:erdc/proteus-docs
  cd proteus
  make docs
If you are  satisified with both the current proteus commit and the docs commit (and push both to github), then tag the release through the github releases site.

- Switch back to dev

Update proteus/__init__.py (add '.dev0' to version)
Update setup.py (add '.dev0' to version)
