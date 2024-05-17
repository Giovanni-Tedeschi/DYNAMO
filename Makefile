dynamo: dynamo.cpp
	g++ -o dynamo.exe dynamo.cpp -lm

clean:
	rm -f dynamo.exe