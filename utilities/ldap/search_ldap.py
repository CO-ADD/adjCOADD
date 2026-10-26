import ldap
#from ldap3 import Server, Connection, ALL, SUBTREE

# 1. Define configuration parameters
LDAP_SERVER = "ldap://ldap.uq.edu.au:636"  # Use ldaps:// for secure connections
# BIND_USER = "CN=Service Account,OU=Users,DC=company,DC=com"
# BIND_PASSWORD = "YourSecurePassword"
BIND_USER = ""
BIND_PASSWORD = ""

SEARCH_BASE = "ou=people,o=The University of Queensland,c=au"
SEARCH_FILTER = "(uid=%(user)s)"
ATTRIBUTES = ["cn", "email", "uid","name"]

# 2. Initialize the server object
# server = Server(LDAP_SERVER, get_info=ALL)

# # 3. Connect and execute search using a context manager
# with Connection(server, user=BIND_USER, password=BIND_PASSWORD, auto_bind=True) as conn:
#     # Perform the search operation
#     status, result, response, _ = conn.search(
#         search_base=SEARCH_BASE,
#         search_filter=SEARCH_FILTER,
#         search_scope=SUBTREE,
#         attributes=ATTRIBUTES
#     )
    
#     # 4. Process the returned entries safely
#     if conn.entries:
#         for entry in conn.entries:
#             print(f"Name: {entry.cn}")
#             print(f"Email: {entry.email}")
#             print(f"Phone: {entry.uid}")
#             print("-" * 20)
#     else:
#         print("No matching entries found.")
        

# Initialize and configure connection
conn = ldap.initialize(LDAP_SERVER)
conn.protocol_version = ldap.VERSION3

# Authenticate/Bind
bind_dn = SEARCH_BASE
password = BIND_PASSWORD
conn.simple_bind_s(bind_dn, password)

# Execute synchronous search
base_dn = SEARCH_BASE
search_filter = "(uid=jdoe)"
attributes = ["mail", "cn"]

try:
    results = conn.search_s(base_dn, ldap.SCOPE_SUBTREE, search_filter, attributes)
    for dn, entry in results:
        if dn: # Ensure it is a valid entry and not a referral
            print(f"DN: {dn}")
            print(f"Attributes: {entry}")
except ldap.LDAPError as e:
    print(f"LDAP Error: {e}")
finally:
    conn.unbind_s()