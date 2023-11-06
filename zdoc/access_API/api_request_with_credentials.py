import requests

token_url = "http://localhost:8000/api/token/"
api_url = "http://localhost:8000/api-drug/"

credentials = {
    "username": "uqName",
    "password": "uqpassword"
}

# Obtain access token
token_response = requests.post(token_url, data=credentials)

if token_response.status_code == 200:
    access_token = token_response.json()["access"]

    # Make an authenticated API request
    headers = {
        "Authorization": f"Bearer {access_token}"
    }

    api_response = requests.get(api_url, headers=headers)

    if api_response.status_code == 200:
        print("Success! Response data:")
        print(api_response.json())
    else:
        print("Error! Status code:", api_response.status_code)
        print("Response data:", api_response.json())
else:
    print("Failed to obtain access token. Status code:", token_response.status_code)
    print("Response data:", token_response.json())