import nose,os
proteus=os.getenv('PROTEUS')
nose.run(defaultTest=os.path.join(proteus,'proteusModule','test'))
print("PEXPECT_EXIT")
