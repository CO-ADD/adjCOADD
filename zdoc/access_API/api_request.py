import requests

api_url = "http://localhost:8000/api-drug/"
token = "generated by python manage.py drf_create_token"

headers = {
    "Authorization": f"Token {token}"
}

response = requests.get(api_url, headers=headers)

if response.status_code == 200:
    print("Success! Response data:")
    # print(response.json())
else:
    print("Error! Status code:", response.status_code)
    print("Response data:", response.json())