rm style_contact_model.whitelist
sh Make.sh models
cp style_contact_model.h style_contact_model.whitelist
gedit style_contact_model.whitelist &

echo "Keep only the entries of models you want to compile."
