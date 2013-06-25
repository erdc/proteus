import nose,os
proteus=os.getenv('PROTEUS')
nose.run(defaultTest=os.path.join(proteus,'externalPackages','petsc4py','test'),argv=['runtests.py','-v','--exclude','test_ts','--exclude','test_da'])
print("PEXPECT_EXIT")
