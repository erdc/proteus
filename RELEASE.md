- To release a new version of proteus
Update proteus/__init__.py (remove '.dev0' and increment version)
Update setup.py (remove 'dev0' and increment version)
Update doc/source/_templates/layout.html (added line for new release in pull down, set to default download to current release)
rebuild docs and check docs
  cd .. 
  git clone https://github.com/erdc-cm/proteus proteus-website -b gh-pages
  cd proteus
  make docs
When you are  satisified with both the current proteus commit and the docs commit and push proteus, and tag the release through the github releases site. Then commit the proteus-website gh-pages branch to publish the documentation on proteustoolkit.org

- Switch back to dev

Update proteus/__init__.py (add '.dev0' to version)
Update setup.py (add '.dev0' to version)
