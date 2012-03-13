import nose,os
proteus=os.getenv('PROTEUS')
nose.run(defaultTest=os.path.join(proteus,'externalPackages','petsc4py','test'),argv=['runtests.py','--exclude','test_ts'])
print("PEXPECT_EXIT")
