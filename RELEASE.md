- To release a new version of proteus

Update doc/source/_templates/layout.html (added line for new release in pull down, set to default download to current release)
rebuild docs and check docs
  cd .. 
  git clone https://github.com/erdc-cm/proteus proteus-website -b gh-pages
  cd proteus
  make doc
When you are  satisified with both the current proteus commit and the docs commit and push proteus, and tag the release through the github releases site. Then commit the proteus-website gh-pages branch to publish the documentation on proteustoolkit.org

- Switch back to dev

Update proteus/__init__.py (add 'dev' to version)
Update setup.py (add 'dev' to version)
