import random
def guess_user_number():
	secret_number = "oo"
	guess_number = 50
	command = ""
	limit_high = 101
	limit_low = 0
	counter = 0
	print("Close your eyes and keep a number in your head between 0-100. No cheating!!!")
	while secret_number != guess_number:
		command = str(input(("Is it {}?(Y) \n(H)igher or (L)ower?\n".format(str(guess_number))))).upper()
		counter += 1
		if command == "H":
			limit_low = guess_number
			guess_number = int((limit_high + guess_number) / 2)
		elif command == "L":
			limit_high = guess_number
			guess_number = int((limit_low + guess_number) / 2)
		elif command == "Y":
			secret_number = guess_number
		else:
			print("Sorry?")
	if counter == 1:
		print("I got you at my first shot! That's how cool I am!")
	else:
		print("I got you in {} tries!!".format(counter))
	command = input("Wanna play again? (Y/N/(Q)uit)").upper()
	if command == "N": return main_program()
	elif command == "Y": return guess_user_number()
	elif command == "Q" or command== "QUIT": return None
	else:
		print("Didn't get you, returning main program")
		return main_program()
def guess_my_number():
	secret_number = random.randint(1,100)
	print("I have a special number for you in my RAM. Can you find it? :)\nHint: it's between 0-100.\nShoot!")
	guess = 0
	guess_count=0
	my_number_command =""
	while guess != secret_number:
		try:
			guess = int(input(""))
		except: TypeError
		if guess > 100 or guess <= 0:
			print("Oh come on! I already gave you the hint!")
		elif guess < secret_number:
			print("Go higher")
		elif guess > secret_number:
			print("Go lower")
		guess_count += 1
	if guess_count == 1:
		print("OMCPU, you just found it in one go!")
	else:
		print("Good you found it! Though it took you {} tries.".format(guess_count))
	my_number_command = input("Wanna play again? (Y/N/(Q)uit)").upper()
	if my_number_command == "N": return main_program()
	elif my_number_command == "Y": return guess_my_number()
	elif my_number_command == "Q" or my_number_command == "QUIT": return None
	else:
		print("Didn't get you, returning main program")
		return main_program()
def main_program():
	maincommand = ""
	print(
		"Welcome to the numbers game! \nM to guess my secret number(Player guesses), Y to guess your secret number(Computer guesses), Q to quit")
	maincommand = str(input("")).upper()
	while maincommand != "Q":
		if maincommand == "M":
			return guess_my_number()
		elif maincommand == "Y":
			return guess_user_number()
		else:
			maincommand = str(input("I didn't quite understand. Again please?")).upper()
	return print("Goodbye!!")

main_program()